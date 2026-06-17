/**
 * @file ntt.h
 * @brief 数论变换（Number Theoretic Transform, NTT）实现
 *
 * NTT 是离散傅里叶变换（DFT）在有限域 Z_q 上的类比。
 * 它是现代格密码方案（如 Kyber/ML-KEM、Dilithium/ML-DSA）的核心算子。
 *
 * 数学定义（Negacyclic NTT）：
 *   给定多项式 a(x) = sum_{i=0}^{n-1} a_i * x^i ∈ Z_q[x] / (x^n + 1)
 *
 *   设 ψ 为 Z_q 中的 2n 次本原单位根（ψ^n ≡ -1 (mod q)），
 *   ω = ψ^2 为 n 次本原单位根（ω^n ≡ 1 (mod q)）。
 *
 *   Negacyclic NTT 正变换:
 *     A_k = sum_{i=0}^{n-1} (a_i * ψ^i) * ω^{ik} mod q
 *         = NTT_standard(a_i * ψ^i)
 *
 *   Negacyclic NTT 逆变换:
 *     a_i = ψ^{-i} * INTT_standard(A_k)_i mod q
 *
 *   性质：negacyclic NTT 域中的逐点相乘对应于
 *   Z_q[x]/(x^n + 1) 上的 negacyclic 卷积（多项式乘法）。
 *
 * 本实现支持两组参数：
 *   1. Kyber 参数:     q = 3329,     n = 256, omega = 17   (n 次单位根)
 *   2. Dilithium 参数: q = 8380417,  n = 256, omega = 1753 (2n 次单位根)
 *
 * 注意：Kyber 的 q-1 = 3328 = 13×2^8，而 2n = 512 = 2^9 不整除 q-1，
 * 因此 Z_q* 中不存在 2n 次单位根 ψ（ψ^n = -1）。预乘/后乘方法不可用。
 * Kyber 实际使用"不完全分裂 NTT"（7 层蝶形 + 二次域 basemul），
 * 本演示工程对 Kyber 的 poly_mul 回退到 O(n²) 朴素乘法以保证正确性。
 * Dilithium 的 q-1 = 8380416 = 2^23×3×7，2n | q-1，可使用完整快速 NTT。
 *
 * @note 本实现侧重教学演示，生产环境请使用高度优化的 NTT 库
 *       （如 libntt, PALISADE/OpenFHE 内置 NTT 等）
 */

#ifndef NTT_H
#define NTT_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

namespace ntt {

// ==========================================================================
// 预定义的密码学参数
// ==========================================================================

/**
 * @brief Kyber/ML-KEM 的 NTT 参数
 *
 * Kyber 使用不完全分裂（incomplete splitting）的 NTT：
 *   q = 3329 (素数，q ≡ 1 mod 256)
 *   n = 256
 *   omega = 17 (Z_q 中的 256 次本原单位根)
 *
 * 验证：17^256 mod 3329 = 1, 17^128 mod 3329 = 3328 ≡ -1 (mod 3329)
 */
struct KyberParams {
    static constexpr uint32_t q     = 3329;     // 模数
    static constexpr uint32_t n     = 256;      // 多项式维度
    static constexpr uint32_t omega = 17;       // 本原单位根
};

/**
 * @brief Dilithium/ML-DSA 的 NTT 参数
 *
 * Dilithium 使用完全分裂（complete splitting）的 NTT：
 *   q = 8380417 (素数，q ≡ 1 mod 512)
 *   n = 256
 *   omega = 1753 (Z_q 中的 512 次本原单位根)
 */
struct DilithiumParams {
    static constexpr uint32_t q     = 8380417;  // 模数
    static constexpr uint32_t n     = 256;      // 多项式维度
    static constexpr uint32_t omega = 1753;     // 本原单位根
};

// ==========================================================================
// 模运算工具函数
// ==========================================================================

/**
 * @brief 安全的模约减，确保结果在 [0, q) 范围内
 * @param a  输入值（可能为负数，使用 int64_t 处理）
 * @param q  模数
 * @return   a mod q，结果在 [0, q) 范围内
 */
inline uint32_t mod_reduce(int64_t a, uint32_t q) {
    int64_t r = a % static_cast<int64_t>(q);
    if (r < 0) r += static_cast<int64_t>(q);
    return static_cast<uint32_t>(r);
}

/**
 * @brief 模加法: (a + b) mod q
 * @param a  操作数 a ∈ [0, q)
 * @param b  操作数 b ∈ [0, q)
 * @param q  模数
 * @return   (a + b) mod q
 */
inline uint32_t mod_add(uint32_t a, uint32_t b, uint32_t q) {
    uint32_t r = a + b;
    if (r >= q) r -= q;
    return r;
}

/**
 * @brief 模减法: (a - b) mod q
 * @param a  操作数 a ∈ [0, q)
 * @param b  操作数 b ∈ [0, q)
 * @param q  模数
 * @return   (a - b + q) mod q
 */
inline uint32_t mod_sub(uint32_t a, uint32_t b, uint32_t q) {
    uint32_t r = a - b;
    if (r > a) r += q;  // 发生了借位
    return r;
}

/**
 * @brief 模乘法: (a * b) mod q
 *
 * 使用 uint64_t 中间结果防止溢出。
 * 对于较大的 q（如 Dilithium 的 q ≈ 2^23），
 * 两个 q 范围内的数相乘最大约 2^46，不会溢出 uint64_t。
 *
 * @param a  操作数 a ∈ [0, q)
 * @param b  操作数 b ∈ [0, q)
 * @param q  模数
 * @return   (a * b) mod q
 */
inline uint32_t mod_mul(uint32_t a, uint32_t b, uint32_t q) {
    return static_cast<uint32_t>(
        (static_cast<uint64_t>(a) * static_cast<uint64_t>(b))
        % static_cast<uint64_t>(q)
    );
}

/**
 * @brief 快速幂（模幂运算）: base^exp mod q
 *
 * 使用二进制展开法（Square-and-Multiply），复杂度 O(log exp)。
 * 这是计算 omega^k mod q 的基础。
 *
 * @param base  底数
 * @param exp   指数
 * @param q     模数
 * @return      base^exp mod q
 */
uint32_t mod_pow(uint32_t base, uint32_t exp, uint32_t q);

/**
 * @brief 模逆元: a^{-1} mod q
 *
 * 使用费马小定理: a^{-1} ≡ a^{q-2} (mod q)（要求 q 为素数）
 * 也可使用扩展欧几里得算法，此处使用快速幂实现。
 *
 * @param a  输入值 ∈ [1, q)
 * @param q  素数模数
 * @return   a^{-1} mod q
 */
uint32_t mod_inv(uint32_t a, uint32_t q);

/**
 * @brief 模平方根: sqrt(a) mod q（Tonelli-Shanks 算法）
 *
 * 计算满足 x^2 ≡ a (mod q) 的 x。
 * 用于从 n 次单位根 omega 计算 2n 次单位根 psi = sqrt(omega)。
 *
 * 当 q ≡ 3 (mod 4) 时可用简化公式 x = a^{(q+1)/4} mod q。
 * 本实现使用通用的 Tonelli-Shanks 算法，适用于所有奇素数 q。
 *
 * @param a  输入值
 * @param q  奇素数模数
 * @return   sqrt(a) mod q，如果不存在则返回 0
 */
uint32_t mod_sqrt(uint32_t a, uint32_t q);

// ==========================================================================
// NTT 核心接口
// ==========================================================================

/**
 * @brief NTT 上下文，预计算蝶形运算所需的旋转因子（twiddle factors）
 *
 * 通过预计算避免在每次 NTT 调用时重复计算 omega 的幂次，
 * 这是实际实现中的常见优化策略。
 *
 * @tparam Params  参数结构体，需包含 q, n, omega 三个静态常量
 */
template <typename Params>
class NTTContext {
public:
    /**
     * @brief 构造 NTT 上下文，预计算 NTT 所需的全部参数
     *
     * 预计算内容：
     *   1. ψ (psi): 2n 次本原单位根（仅当存在时）
     *      - 若 omega^n ≡ -1 (mod q): ψ = omega（Dilithium，快速 NTT 可用）
     *      - 若 omega^n ≡ 1 且 2n 不整除 q-1: ψ 不存在（Kyber，回退朴素乘法）
     *   2. ω_ntt = ψ^2: 标准 NTT 使用的 n 次单位根
     *   3. 正变换旋转因子: ω_ntt^i mod q
     *   4. 逆变换旋转因子: ω_ntt^{-i} mod q
     *   5. ψ^{-1} 和 n^{-1}（用于逆变换的归一化）
     */
    NTTContext();

    /**
     * @brief NTT 正向变换（Cooley-Tukey 蝶形）
     *
     * 当 use_fast_ntt_ = true（如 Dilithium）时执行完整 negacyclic NTT：
     *   1. 预乘: poly[i] = poly[i] * ψ^i mod q（negacyclic twist）
     *   2. 标准 NTT: 比特反转 + Cooley-Tukey 蝶形（使用 ω_ntt = ψ^2）
     *
     * 当 use_fast_ntt_ = false（如 Kyber）时仅执行标准 NTT（无预乘）。
     */
    void forward(uint32_t* poly) const;

    /**
     * @brief NTT 逆变换（Gentleman-Sande 蝶形）
     *
     * 当 use_fast_ntt_ = true（如 Dilithium）时执行完整 negacyclic INTT：
     *   1. 标准 INTT: GS 蝶形 + 比特反转 + 归一化（使用 ω_ntt^{-1}）
     *   2. 后乘: poly[i] = poly[i] * ψ^{-i} mod q（撤销 negacyclic twist）
     *
     * 当 use_fast_ntt_ = false（如 Kyber）时仅执行标准 INTT（无后乘）。
     */
    void inverse(uint32_t* poly) const;

    /**
     * @brief NTT 域多项式乘法: c = a ⊙ b (逐点相乘)
     *
     * 在 NTT 域中，多项式乘法简化为逐点相乘：
     *   c_i = a_i * b_i mod q,  i = 0,...,n-1
     *
     * 这是 NTT 的核心价值——将卷积运算转化为逐点运算。
     *
     * @param a  NTT 域多项式 a（不被修改）
     * @param b  NTT 域多项式 b（不被修改）
     * @param c  输出多项式 c（NTT 域）
     */
    void pointwise_mul(const uint32_t* a, const uint32_t* b, uint32_t* c) const;

    /**
     * @brief 完整的多项式乘法流程: c = a * b in Z_q[x]/(x^n+1)
     *
     * 步骤：
     *   1. a_ntt = NTT(a)
     *   2. b_ntt = NTT(b)
     *   3. c_ntt = a_ntt ⊙ b_ntt（逐点相乘）
     *   4. c = INTT(c_ntt)
     *
     * 总复杂度: O(n log n)，远优于朴素乘法的 O(n^2)
     *
     * @param a  输入多项式 a（系数表示，不被修改）
     * @param b  输入多项式 b（系数表示，不被修改）
     * @param c  输出多项式 c（系数表示）
     */
    void poly_mul(const uint32_t* a, const uint32_t* b, uint32_t* c) const;

    /**
     * @brief 获取模数 q
     */
    uint32_t get_modulus() const { return q_; }

    /**
     * @brief 获取多项式维度 n
     */
    uint32_t get_degree() const { return n_; }

private:
    uint32_t q_;                        // 模数
    uint32_t n_;                        // 多项式维度
    uint32_t psi_;                      // 2n 次本原单位根（ψ^n ≡ -1 mod q）
    uint32_t psi_inv_;                  // ψ 的模逆元（用于 INTT 后乘）
    uint32_t omega_ntt_;                // NTT 使用的 n 次单位根（= ψ^2）
    uint32_t omega_ntt_inv_;            // ω_ntt 的模逆元
    uint32_t n_inv_;                    // n 的模逆元（用于 INTT 归一化）
    bool use_fast_ntt_;                 // 是否可使用快速 NTT（需要 ψ 存在）
    std::vector<uint32_t> twiddles_;    // 正变换旋转因子: ω_ntt^i
    std::vector<uint32_t> inv_twiddles_;// 逆变换旋转因子: ω_ntt^{-i}

    /**
     * @brief 比特反转置换（Bit-Reversal Permutation）
     *
     * NTT 的迭代实现需要先将输入按比特反转顺序排列。
     * 例如 n=8 时：
     *   索引 0(000) -> 0(000)
     *   索引 1(001) -> 4(100)
     *   索引 2(010) -> 2(010)
     *   索引 3(011) -> 6(110)
     *   ...
     *
     * @param poly  输入/输出数组（原地置换）
     */
    void bit_reverse(uint32_t* poly) const;

    /**
     * @brief 计算 log2(n)（n 必须是 2 的幂）
     */
    uint32_t log2n() const;
};

// ==========================================================================
// 辅助函数
// ==========================================================================

/**
 * @brief 验证 NTT 参数的正确性
 *
 * 检查：
 *   1. q 是否为素数（简单的试除法）
 *   2. n 是否为 2 的幂
 *   3. omega 是否为本原单位根
 *      - omega^n ≡ -1: 2n 次单位根（Dilithium，快速 NTT 可用）
 *      - omega^n ≡ 1 且 omega^{n/2} ≡ -1: n 次单位根（Kyber，回退朴素乘法）
 *
 * @tparam Params  参数结构体
 * @return true 表示参数正确，false 表示参数有误
 */
template <typename Params>
bool verify_params();

/**
 * @brief 打印多项式系数（用于调试）
 * @param name  多项式名称
 * @param poly  系数数组
 * @param n     多项式维度
 * @param q     模数
 */
void print_poly(const std::string& name, const uint32_t* poly,
                size_t n, uint32_t q);

/**
 * @brief 朴素多项式乘法（用于与 NTT 结果对比验证）
 *
 * 在 Z_q[x]/(x^n + 1) 上执行 O(n^2) 的朴素乘法。
 * 仅用于小规模测试和正确性验证。
 *
 * @param a  多项式 a
 * @param b  多项式 b
 * @param c  输出多项式 c = a * b mod (x^n + 1, q)
 * @param n  多项式维度
 * @param q  模数
 */
void naive_poly_mul(const uint32_t* a, const uint32_t* b, uint32_t* c,
                    size_t n, uint32_t q);

} // namespace ntt

#endif // NTT_H
