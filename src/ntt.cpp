/**
 * @file ntt.cpp
 * @brief 数论变换（NTT）实现 — 面向格密码的教学级实现
 *
 * 本文件实现了完整的 NTT 正/逆变换及其辅助函数。
 * 实现采用迭代 Cooley-Tukey / Gentleman-Sande 蝶形算法，
 * 避免了递归调用的开销，同时保持代码的可读性。
 *
 * 参考文献：
 *   [1] Lyubashevsky et al., "CRYSTALS-Kyber", NIST PQC Standard
 *   [2] Ducas et al., "CRYSTALS-Dilithium", NIST PQC Standard
 *   [3] Seiler, "NTT: Number Theoretic Transform", 2023
 *   [4] Harvey, "Faster truncated polynomial multiplication", 2019
 */

#include "ntt.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>

namespace ntt {

// ==========================================================================
// 模运算工具函数实现
// ==========================================================================

uint32_t mod_pow(uint32_t base, uint32_t exp, uint32_t q)
{
    /**
     * 快速幂算法（Binary Exponentiation / Square-and-Multiply）
     *
     * 原理：将指数 exp 按二进制展开
     *   base^exp = base^{b_k * 2^k + ... + b_1 * 2 + b_0}
     *
     * 从高位到低位扫描：
     *   result = 1
     *   for each bit from MSB to LSB:
     *     result = result^2 mod q         （平方）
     *     if bit == 1:
     *       result = result * base mod q  （乘入）
     *
     * 复杂度：O(log exp) 次模乘法
     *
     * 注意：此处使用 uint64_t 存储中间结果以防止溢出。
     *       当 q < 2^31 时，两个 < q 的数相乘 < 2^62，不会溢出 uint64_t。
     */
    uint64_t result = 1;
    uint64_t b = static_cast<uint64_t>(base) % static_cast<uint64_t>(q);

    while (exp > 0) {
        // 如果当前位为 1，乘入 base
        if (exp & 1) {
            result = (result * b) % static_cast<uint64_t>(q);
        }
        // 平方 base
        b = (b * b) % static_cast<uint64_t>(q);
        exp >>= 1;
    }

    return static_cast<uint32_t>(result);
}

uint32_t mod_inv(uint32_t a, uint32_t q)
{
    /**
     * 模逆元计算（基于费马小定理）
     *
     * 费马小定理：当 q 为素数且 gcd(a, q) = 1 时，
     *   a^{q-1} ≡ 1 (mod q)
     *   => a^{-1} ≡ a^{q-2} (mod q)
     *
     * 替代方法（性能更优，此处未采用）：
     *   - 扩展欧几里得算法（Extended Euclidean Algorithm）
     *   - 对于固定的 q，可以预计算 Montgomery 参数
     */
    assert(a > 0 && "模逆元的输入不能为零");
    return mod_pow(a, q - 2, q);
}

uint32_t mod_sqrt(uint32_t a, uint32_t q)
{
    /**
     * Tonelli-Shanks 算法：计算 sqrt(a) mod q
     *
     * 算法步骤：
     *   1. 分解 q - 1 = Q * 2^S（提取 2 的幂次）
     *   2. 找到一个模 q 的二次非剩余 z
     *   3. 初始化：M = S, c = z^Q, t = a^Q, R = a^{(Q+1)/2}
     *   4. 循环：找到最小的 i 使得 t^{2^i} = 1
     *      - 若 i = 0，则 R 即为所求
     *      - 否则更新 M, c, t, R 并继续
     *
     * 复杂度：O(log^2 q)
     *
     * 在 NTT 中的应用：
     *   当 omega 是 n 次单位根（omega^n = 1）时，
     *   需要计算 psi = sqrt(omega) 作为 2n 次单位根，
     *   用于 negacyclic NTT 的预乘/后乘步骤。
     */

    if (a == 0) return 0;
    if (q == 2) return a;

    // 检查 a 是否为二次剩余（欧拉准则）
    if (mod_pow(a, (q - 1) / 2, q) != 1) {
        return 0;  // a 不是模 q 的二次剩余
    }

    // 分解 q - 1 = Q * 2^S
    uint32_t Q = q - 1;
    uint32_t S = 0;
    while (Q % 2 == 0) {
        Q /= 2;
        S++;
    }

    // 特殊情况：S = 1（即 q ≡ 3 mod 4），可用简化公式
    if (S == 1) {
        return mod_pow(a, (q + 1) / 4, q);
    }

    // 找到二次非剩余 z
    uint32_t z = 2;
    while (mod_pow(z, (q - 1) / 2, q) != q - 1) {
        z++;
    }

    uint32_t M = S;
    uint32_t c = mod_pow(z, Q, q);
    uint32_t t = mod_pow(a, Q, q);
    uint32_t R = mod_pow(a, (Q + 1) / 2, q);

    while (true) {
        if (t == 1) return R;

        // 找到最小的 i (0 < i < M) 使得 t^{2^i} ≡ 1 (mod q)
        uint32_t i = 0;
        uint32_t tmp = t;
        do {
            tmp = mod_mul(tmp, tmp, q);
            i++;
        } while (tmp != 1 && i < M);

        // 更新参数
        uint32_t b = c;
        for (uint32_t j = 0; j < M - i - 1; ++j) {
            b = mod_mul(b, b, q);
        }
        M = i;
        c = mod_mul(b, b, q);
        t = mod_mul(t, c, q);
        R = mod_mul(R, b, q);
    }
}

// ==========================================================================
// NTTContext 模板实现
// ==========================================================================

template <typename Params>
NTTContext<Params>::NTTContext()
    : q_(Params::q)
    , n_(Params::n)
{
    /**
     * 构造函数：推导 negacyclic NTT 所需的全部参数
     *
     * 核心步骤：
     *   1. 从模板参数 omega 推导 ψ（2n 次本原单位根）
     *   2. 计算 ω_ntt = ψ^2（标准 NTT 的 n 次单位根）
     *   3. 预计算旋转因子和逆旋转因子
     *
     * 两种参数情况的处理：
     *
     *   情况 A（Dilithium）：omega^n ≡ -1 (mod q)
     *     omega 本身即为 2n 次单位根，直接令 ψ = omega
     *     ω_ntt = omega^2（n 次单位根）
     *
     *   情况 B（Kyber）：omega^n ≡ 1 (mod q)
     *     omega 是 n 次单位根，需要计算 ψ = sqrt(omega)
     *     使用 Tonelli-Shanks 算法求模平方根
     *     ω_ntt = omega（n 次单位根不变）
     */

    constexpr uint32_t omega = Params::omega;

    // 计算 n 的逆元（用于 INTT 的归一化）
    n_inv_ = mod_inv(n_, q_);

    // 判断 omega 是 n 次还是 2n 次单位根
    uint32_t omega_n = mod_pow(omega, n_, q_);

    if (omega_n == q_ - 1) {
        // 情况 A：omega^n = -1，omega 本身是 2n 次单位根（Dilithium）
        psi_ = omega;
        omega_ntt_ = mod_mul(omega, omega, q_);  // ω_ntt = ψ^2
    } else {
        // 情况 B：omega^n = 1，omega 是 n 次单位根，需计算 ψ = sqrt(omega)（Kyber）
        psi_ = mod_sqrt(omega, q_);
        assert(psi_ != 0 && "sqrt(omega) 不存在，参数有误");
        // 验证 ψ^n ≡ -1 (mod q)
        assert(mod_pow(psi_, n_, q_) == q_ - 1 && "ψ^n ≠ -1 (mod q)，参数有误");
        omega_ntt_ = omega;  // ω_ntt 不变
    }

    // 计算 ψ 和 ω_ntt 的逆元
    psi_inv_ = mod_inv(psi_, q_);
    omega_ntt_inv_ = mod_inv(omega_ntt_, q_);

    // 分配并预计算旋转因子
    twiddles_.resize(n_);
    inv_twiddles_.resize(n_);

    uint32_t w_pow = 1;       // ω_ntt^0
    uint32_t w_inv_pow = 1;   // ω_ntt^{-0}

    for (uint32_t i = 0; i < n_; ++i) {
        twiddles_[i] = w_pow;
        inv_twiddles_[i] = w_inv_pow;
        w_pow = mod_mul(w_pow, omega_ntt_, q_);
        w_inv_pow = mod_mul(w_inv_pow, omega_ntt_inv_, q_);
    }
}

template <typename Params>
uint32_t NTTContext<Params>::log2n() const
{
    /**
     * 计算 log2(n)，n 必须是 2 的幂
     * 例如 n=256 -> log2n=8
     */
    uint32_t log = 0;
    uint32_t tmp = n_;
    while (tmp > 1) {
        tmp >>= 1;
        log++;
    }
    return log;
}

template <typename Params>
void NTTContext<Params>::bit_reverse(uint32_t* poly) const
{
    /**
     * 比特反转置换（Bit-Reversal Permutation）
     *
     * 迭代 NTT 需要输入数据按比特反转顺序排列。
     * 例如 n=8 (log2n=3) 时：
     *
     *   原始索引  二进制  反转后  目标索引
     *      0      000     000      0
     *      1      001     100      4
     *      2      010     010      2
     *      3      011     110      6
     *      4      100     001      1
     *      5      101     101      5
     *      6      110     011      3
     *      7      111     111      7
     *
     * 实现优化：只交换 i < rev(i) 的对，避免重复交换
     */
    uint32_t log = log2n();

    for (uint32_t i = 0; i < n_; ++i) {
        // 计算 i 的比特反转
        uint32_t rev = 0;
        uint32_t tmp = i;
        for (uint32_t bit = 0; bit < log; ++bit) {
            rev = (rev << 1) | (tmp & 1);
            tmp >>= 1;
        }

        // 仅当 i < rev 时交换，避免重复交换
        if (i < rev) {
            std::swap(poly[i], poly[rev]);
        }
    }
}

template <typename Params>
void NTTContext<Params>::forward(uint32_t* poly) const
{
    /**
     * Negacyclic NTT 正向变换
     *
     * 步骤：
     *   1. 预乘 ψ^i（negacyclic twist）：poly[i] *= ψ^i mod q
     *   2. 标准 NTT（Cooley-Tukey 蝶形）：比特反转 + 蝶形运算
     *
     * 预乘 ψ^i 的数学意义：
     *   将 Z_q[x]/(x^n+1) 上的 negacyclic 卷积映射到
     *   标准 NTT 可处理的 cyclic 卷积。
     *   ψ 是 2n 次单位根（ψ^n = -1），确保了 x^n ≡ -1 的约减规则。
     *
     * 蝶形操作（标准 CT）：
     *   t = ω_ntt^{j*step} * poly[k + j + half] mod q
     *   poly[k + j + half] = poly[k + j] - t mod q
     *   poly[k + j]        = poly[k + j] + t mod q
     */

    // 第一步：预乘 ψ^i（negacyclic twist）
    uint32_t psi_pow = 1;  // ψ^0 = 1
    for (uint32_t i = 0; i < n_; ++i) {
        poly[i] = mod_mul(poly[i], psi_pow, q_);
        psi_pow = mod_mul(psi_pow, psi_, q_);
    }

    // 第二步：比特反转置换
    bit_reverse(poly);

    uint32_t log = log2n();

    // 逐轮蝶形运算
    for (uint32_t s = 1; s <= log; ++s) {
        uint32_t len = 1u << s;           // 当前子块长度: 2, 4, 8, ..., n
        uint32_t half = len >> 1;          // 子块的一半: 1, 2, 4, ..., n/2
        uint32_t step = n_ / len;          // 旋转因子的步长

        // 遍历所有子块
        for (uint32_t k = 0; k < n_; k += len) {
            // 遍历子块内的蝶形对
            for (uint32_t j = 0; j < half; ++j) {
                // 获取旋转因子 w = omega^{j * step}
                uint32_t w = twiddles_[j * step];

                // 蝶形操作
                uint32_t u = poly[k + j];                     // 上支路
                uint32_t t = mod_mul(w, poly[k + j + half], q_); // 下支路 × 旋转因子

                poly[k + j]        = mod_add(u, t, q_);       // 加法
                poly[k + j + half] = mod_sub(u, t, q_);       // 减法
            }
        }
    }
}

template <typename Params>
void NTTContext<Params>::inverse(uint32_t* poly) const
{
    /**
     * Negacyclic NTT 逆向变换
     *
     * 步骤：
     *   1. 标准 INTT（Gentleman-Sande 蝶形）：逆蝶形 + 比特反转 + 归一化
     *   2. 后乘 ψ^{-i}（撤销 negacyclic twist）：poly[i] *= ψ^{-i} mod q
     *
     * GS 蝶形（单步）：
     *   u = poly[k + j]
     *   v = poly[k + j + half]
     *   poly[k + j]        = u + v mod q
     *   poly[k + j + half] = ω_ntt^{-j*step} * (u - v) mod q
     */

    uint32_t log = log2n();

    // 逐轮蝶形运算（方向与正向相反）
    for (uint32_t s = log; s >= 1; --s) {
        uint32_t len = 1u << s;
        uint32_t half = len >> 1;
        uint32_t step = n_ / len;

        for (uint32_t k = 0; k < n_; k += len) {
            for (uint32_t j = 0; j < half; ++j) {
                // 使用逆变换旋转因子
                uint32_t w = inv_twiddles_[j * step];

                uint32_t u = poly[k + j];
                uint32_t v = poly[k + j + half];

                // GS 蝶形：先加减，再乘旋转因子
                poly[k + j]        = mod_add(u, v, q_);
                poly[k + j + half] = mod_mul(w, mod_sub(u, v, q_), q_);
            }
        }
    }

    // 比特反转置换（GS 蝶形在逆序输入上操作，输出为自然序）
    bit_reverse(poly);

    // 归一化：每个系数乘以 n^{-1} mod q
    for (uint32_t i = 0; i < n_; ++i) {
        poly[i] = mod_mul(poly[i], n_inv_, q_);
    }

    // 后乘 ψ^{-i}（撤销 negacyclic twist）
    // 正变换预乘了 ψ^i，逆变换需要乘以 ψ^{-i} 来恢复原始系数
    uint32_t psi_inv_pow = 1;  // ψ^{-0} = 1
    for (uint32_t i = 0; i < n_; ++i) {
        poly[i] = mod_mul(poly[i], psi_inv_pow, q_);
        psi_inv_pow = mod_mul(psi_inv_pow, psi_inv_, q_);
    }
}

template <typename Params>
void NTTContext<Params>::pointwise_mul(const uint32_t* a, const uint32_t* b,
                                        uint32_t* c) const
{
    /**
     * NTT 域逐点乘法
     *
     * c[i] = a[i] * b[i] mod q,  i = 0, ..., n-1
     *
     * 这是 NTT 的核心价值所在：
     *   多项式乘法（卷积）在 NTT 域中变为逐点相乘
     *   时间复杂度从 O(n^2) 降至 O(n)
     *
     * 在实际的格密码实现中，这一步通常是最快的，
     * 性能瓶颈往往在 NTT/INTT 变换本身。
     */
    for (uint32_t i = 0; i < n_; ++i) {
        c[i] = mod_mul(a[i], b[i], q_);
    }
}

template <typename Params>
void NTTContext<Params>::poly_mul(const uint32_t* a, const uint32_t* b,
                                   uint32_t* c) const
{
    /**
     * 完整的多项式乘法流程（使用 NTT 加速）
     *
     * 步骤：
     *   1. 复制输入（避免修改原始数据）
     *   2. a_ntt = NTT(a)     —— 正向变换
     *   3. b_ntt = NTT(b)     —— 正向变换
     *   4. c_ntt = a_ntt ⊙ b_ntt  —— 逐点相乘
     *   5. c = INTT(c_ntt)    —— 逆向变换
     *
     * 总复杂度: 2 * O(n log n) + O(n) = O(n log n)
     * 对比朴素乘法: O(n^2)
     *
     * 对于 Kyber 参数 (n=256):
     *   NTT 加速比 ≈ n / (2 log n) ≈ 256 / 16 = 16 倍
     *
     * 对于 Dilithium 参数 (n=256, 但 q 更大):
     *   加速比类似，但模运算开销更大
     */

    // 分配临时数组
    std::vector<uint32_t> a_ntt(n_);
    std::vector<uint32_t> b_ntt(n_);

    // 复制输入
    std::memcpy(a_ntt.data(), a, n_ * sizeof(uint32_t));
    std::memcpy(b_ntt.data(), b, n_ * sizeof(uint32_t));

    // 正向 NTT
    forward(a_ntt.data());
    forward(b_ntt.data());

    // 逐点相乘
    pointwise_mul(a_ntt.data(), b_ntt.data(), c);

    // 逆向 NTT
    inverse(c);
}

// ==========================================================================
// 辅助函数实现
// ==========================================================================

template <typename Params>
bool verify_params()
{
    /**
     * 验证 NTT 参数的正确性（Negacyclic NTT 版本）
     *
     * 检查项：
     *   1. n 是否为 2 的幂
     *   2. q 是否为素数（简单试除法）
     *   3. 从 omega 推导 ψ（2n 次本原单位根），验证 ψ^n ≡ -1 (mod q)
     *      - omega^n ≡ -1: ψ = omega（Dilithium 情况）
     *      - omega^n ≡ 1:  ψ = sqrt(omega)（Kyber 情况）
     */
    constexpr uint32_t q = Params::q;
    constexpr uint32_t n = Params::n;
    constexpr uint32_t omega = Params::omega;

    // 检查 n 是 2 的幂
    if ((n & (n - 1)) != 0) {
        std::cerr << "错误: n=" << n << " 不是 2 的幂\n";
        return false;
    }

    // 简单的素性检验（试除法）
    if (q < 2) {
        std::cerr << "错误: q=" << q << " 不是素数\n";
        return false;
    }
    for (uint32_t i = 2; static_cast<uint64_t>(i) * i <= q; ++i) {
        if (q % i == 0) {
            std::cerr << "错误: q=" << q << " 不是素数（被 " << i << " 整除）\n";
            return false;
        }
    }

    // 推导 ψ 并验证 ψ^n ≡ -1 (mod q)
    uint32_t psi;
    uint32_t omega_n = mod_pow(omega, n, q);

    if (omega_n == q - 1) {
        // omega 本身是 2n 次单位根（Dilithium 情况）
        psi = omega;
    } else if (omega_n == 1) {
        // omega 是 n 次单位根，需计算 ψ = sqrt(omega)（Kyber 情况）
        psi = mod_sqrt(omega, q);
        if (psi == 0) {
            std::cerr << "错误: sqrt(" << omega << ") 在 Z_" << q << " 中不存在\n";
            return false;
        }
    } else {
        std::cerr << "错误: omega^n = " << omega_n
                  << "，既不是 1 也不是 -1 (mod " << q << ")\n";
        return false;
    }

    // 最终验证：ψ^n ≡ -1 (mod q)
    uint32_t psi_n = mod_pow(psi, n, q);
    if (psi_n != q - 1) {
        std::cerr << "错误: ψ^n = " << psi_n << " != " << (q - 1)
                  << " (mod " << q << ")\n";
        return false;
    }

    return true;
}

void print_poly(const std::string& name, const uint32_t* poly,
                size_t n, uint32_t q)
{
    std::cout << name << " (mod " << q << ", deg < " << n << "):\n  ";
    for (size_t i = 0; i < n; ++i) {
        std::cout << poly[i];
        if (i < n - 1) std::cout << ", ";
        // 每 16 个系数换行，便于阅读
        if ((i + 1) % 16 == 0 && i < n - 1) {
            std::cout << "\n  ";
        }
    }
    std::cout << "\n" << std::endl;
}

void naive_poly_mul(const uint32_t* a, const uint32_t* b, uint32_t* c,
                    size_t n, uint32_t q)
{
    /**
     * 朴素多项式乘法: c = a * b mod (x^n + 1, q)
     *
     * 在环 Z_q[x]/(x^n + 1) 上：
     *   x^n ≡ -1 (mod x^n + 1)
     *
     * 因此，当乘积的次数 >= n 时，需要约减：
     *   x^{n+k} = x^k * x^n ≡ -x^k (mod x^n + 1)
     *
     * 算法复杂度: O(n^2)
     * 仅用于验证 NTT 结果的正确性
     */
    // 清零输出
    std::memset(c, 0, n * sizeof(uint32_t));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // a[i] * b[j] * x^{i+j}
            uint32_t coeff = mod_mul(a[i], b[j], q);
            size_t deg = i + j;

            if (deg < n) {
                // 次数 < n，直接累加
                c[deg] = mod_add(c[deg], coeff, q);
            } else {
                // 次数 >= n，利用 x^n ≡ -1 约减
                // x^{deg} = x^{deg-n} * x^n ≡ -x^{deg-n}
                c[deg - n] = mod_sub(c[deg - n], coeff, q);
            }
        }
    }
}

// ==========================================================================
// 模板显式实例化
// ==========================================================================
// 为 Kyber 和 Dilithium 参数生成具体的模板实例

template class NTTContext<KyberParams>;
template class NTTContext<DilithiumParams>;

template bool verify_params<KyberParams>();
template bool verify_params<DilithiumParams>();

} // namespace ntt
