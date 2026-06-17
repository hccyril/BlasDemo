/**
 * @file main.cpp
 * @brief OpenBLAS 演示项目 — 主程序入口
 *
 * 本程序演示以下核心功能：
 *   1. BLAS Level 1: 向量操作（点积、AXPY、范数等）
 *   2. BLAS Level 2: 矩阵-向量乘法 (GEMV)
 *   3. BLAS Level 3: 矩阵乘法 (GEMM) — OpenBLAS 的核心优势
 *   4. NTT: 数论变换及其在多项式乘法中的应用
 *
 * 运行方式:
 *   ./openblas_demo              # 运行全部演示
 *   ./openblas_demo --blas       # 仅运行 BLAS 演示
 *   ./openblas_demo --ntt        # 仅运行 NTT 演示
 *   ./openblas_demo --bench      # 运行矩阵乘法性能基准测试
 *
 * 编译方法详见 README.md
 */

#include "blas_demo.h"
#include "ntt.h"

// OpenBLAS 头文件（用于获取版本信息和线程控制）
extern "C" {
#include <cblas.h>
#include <openblas_config.h>
}

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

// ==========================================================================
// 辅助宏与函数
// ==========================================================================

/**
 * @brief 计时工具 — 使用高精度时钟测量代码执行时间
 */
class Timer {
public:
    void start() { start_ = std::chrono::high_resolution_clock::now(); }

    double elapsed_ms() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }

    double elapsed_us() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::micro>(end - start_).count();
    }

private:
    std::chrono::high_resolution_clock::time_point start_;
};

/// 打印分隔线
void print_separator(const std::string& title)
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n\n";
}

// ==========================================================================
// BLAS Level 1 演示
// ==========================================================================

void demo_blas_level1()
{
    print_separator("BLAS Level 1: 向量-向量操作");

    const size_t n = 8;
    std::vector<double> x(n), y(n);

    // 初始化向量
    blas_demo::random_vector(x.data(), n, -2.0, 2.0);
    blas_demo::random_vector(y.data(), n, -2.0, 2.0);

    blas_demo::print_vector("x", x.data(), n);
    blas_demo::print_vector("y", y.data(), n);

    // --- 1. 向量点积 ---
    std::cout << "--- 1. 向量点积 (DDOT) ---\n";
    double dot = blas_demo::dot_product(x.data(), y.data(), n);
    std::cout << "x^T * y = " << std::fixed << std::setprecision(6) << dot << "\n\n";

    // --- 2. AXPY 操作 ---
    std::cout << "--- 2. AXPY 操作: y = 2.0 * x + y ---\n";
    std::vector<double> y_copy(y);
    blas_demo::axpy(2.0, x.data(), y_copy.data(), n);
    blas_demo::print_vector("y (after axpy)", y_copy.data(), n);

    // --- 3. 向量缩放 ---
    std::cout << "--- 3. 向量缩放 (SCAL): x = 0.5 * x ---\n";
    std::vector<double> x_copy(x);
    blas_demo::scal(0.5, x_copy.data(), n);
    blas_demo::print_vector("x (after scal)", x_copy.data(), n);

    // --- 4. 向量范数 ---
    std::cout << "--- 4. 向量范数 ---\n";
    double l2 = blas_demo::vector_norm(x.data(), n);
    double l1 = blas_demo::vector_asum(x.data(), n);
    std::cout << "||x||_2 (L2 范数) = " << l2 << "\n";
    std::cout << "||x||_1 (L1 范数) = " << l1 << "\n\n";
}

// ==========================================================================
// BLAS Level 2 演示
// ==========================================================================

void demo_blas_level2()
{
    print_separator("BLAS Level 2: 矩阵-向量乘法 (GEMV)");

    const size_t m = 4, n = 3;
    std::vector<double> A(m * n);  // m x n 矩阵（列优先存储）
    std::vector<double> x(n);
    std::vector<double> y(m);

    // 初始化矩阵和向量
    blas_demo::random_matrix(A.data(), m, n, -1.0, 1.0);
    blas_demo::random_vector(x.data(), n, -2.0, 2.0);

    blas_demo::print_matrix("A", A.data(), m, n, m);
    blas_demo::print_vector("x", x.data(), n);

    // 计算 y = A * x
    std::cout << "--- 计算 y = A * x ---\n";
    blas_demo::gemv(A.data(), x.data(), y.data(), m, n, m);
    blas_demo::print_vector("y = A * x", y.data(), m);

    // 计算 y = 2.0 * A * x + 1.0 * y（alpha=2.0, beta=1.0）
    std::cout << "--- 计算 y = 2.0 * A * x + 1.0 * y ---\n";
    blas_demo::gemv(A.data(), x.data(), y.data(), m, n, m, 2.0, 1.0);
    blas_demo::print_vector("y (updated)", y.data(), m);
}

// ==========================================================================
// BLAS Level 3 演示
// ==========================================================================

void demo_blas_level3()
{
    print_separator("BLAS Level 3: 矩阵乘法 (GEMM) — OpenBLAS 核心");

    const size_t m = 4, n = 3, k = 5;
    std::vector<double> A(m * k);  // m x k
    std::vector<double> B(k * n);  // k x n
    std::vector<double> C(m * n);  // m x n

    // 初始化矩阵
    blas_demo::random_matrix(A.data(), m, k, -1.0, 1.0);
    blas_demo::random_matrix(B.data(), k, n, -1.0, 1.0);

    blas_demo::print_matrix("A", A.data(), m, k, m);
    blas_demo::print_matrix("B", B.data(), k, n, k);

    // --- 1. 基本矩阵乘法: C = A * B ---
    std::cout << "--- 1. 计算 C = A * B ---\n";
    blas_demo::gemm(A.data(), B.data(), C.data(), m, n, k, m, k, m);
    blas_demo::print_matrix("C = A * B", C.data(), m, n, m);

    // --- 2. 带缩放的矩阵乘法: C = 2.0 * A * B + 0.5 * C ---
    std::cout << "--- 2. 计算 C = 2.0 * A * B + 0.5 * C ---\n";
    blas_demo::gemm(A.data(), B.data(), C.data(), m, n, k, m, k, m, 2.0, 0.5);
    blas_demo::print_matrix("C (updated)", C.data(), m, n, m);

    // --- 3. 对称矩阵乘法 ---
    std::cout << "--- 3. 对称矩阵乘法 (SYMM) ---\n";
    const size_t sm = 3;
    std::vector<double> S(sm * sm);  // 3x3 对称矩阵
    std::vector<double> B2(sm * 2);  // 3x2 矩阵
    std::vector<double> C2(sm * 2);  // 3x2 输出

    // 构造对称矩阵：先随机填充，再对称化
    blas_demo::random_matrix(S.data(), sm, sm, -1.0, 1.0);
    // 使 S 对称：S[j*sm+i] = S[i*sm+j]（列优先）
    for (size_t i = 0; i < sm; ++i) {
        for (size_t j = i + 1; j < sm; ++j) {
            S[j * sm + i] = S[i * sm + j];
        }
    }
    blas_demo::random_matrix(B2.data(), sm, 2, -1.0, 1.0);

    blas_demo::print_matrix("S (对称)", S.data(), sm, sm, sm);
    blas_demo::print_matrix("B2", B2.data(), sm, 2, sm);

    blas_demo::symm(S.data(), B2.data(), C2.data(), sm, 2, sm, sm, sm, true);
    blas_demo::print_matrix("C2 = S * B2", C2.data(), sm, 2, sm);
}

// ==========================================================================
// GEMM 性能基准测试
// ==========================================================================

void benchmark_gemm()
{
    print_separator("性能基准测试: DGEMM 矩阵乘法（OpenBLAS vs 朴素实现）");

    // 测试不同规模的矩阵
    std::vector<size_t> sizes = {128, 256, 512, 1024, 2048};

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "矩阵规模      | OpenBLAS (ms) | 朴素实现 (ms) | 加速比  | 一致\n";
    std::cout << std::string(72, '-') << "\n";

    for (size_t n : sizes) {
        std::vector<double> A(n * n), B(n * n), C_blas(n * n), C_naive(n * n);
        blas_demo::random_matrix(A.data(), n, n);
        blas_demo::random_matrix(B.data(), n, n);

        // === OpenBLAS GEMM ===
        // 预热（第一次调用可能有内存分配开销）
        blas_demo::gemm(A.data(), B.data(), C_blas.data(), n, n, n, n, n, n);

        const int repeat = (n <= 512) ? 10 : 3;
        Timer timer;
        timer.start();
        for (int i = 0; i < repeat; ++i) {
            blas_demo::gemm(A.data(), B.data(), C_blas.data(), n, n, n, n, n, n);
        }
        double blas_time = timer.elapsed_ms() / repeat;

        // === 朴素 GEMM ===
        // 预热
        blas_demo::naive_gemm(A.data(), B.data(), C_naive.data(), n, n, n, n, n, n);

        // 小规模矩阵多跑几次取平均，大规模只跑一次
        const int naive_repeat = (n <= 256) ? 3 : 1;
        timer.start();
        for (int i = 0; i < naive_repeat; ++i) {
            blas_demo::naive_gemm(A.data(), B.data(), C_naive.data(), n, n, n, n, n, n);
        }
        double naive_time = timer.elapsed_ms() / naive_repeat;

        // 计算加速比
        double speedup = (blas_time > 0.001) ? naive_time / blas_time : 0.0;

        // 比较两种实现的结果是否一致
        // 浮点运算的累加顺序不同会导致微小差异，使用相对容差判断
        bool match = true;
        const double eps = 1e-8;  // 相对容差
        for (size_t idx = 0; idx < n * n; ++idx) {
            double diff = std::abs(C_blas[idx] - C_naive[idx]);
            double scale = std::max({std::abs(C_blas[idx]), std::abs(C_naive[idx]), 1.0});
            if (diff / scale > eps) {
                match = false;
                break;
            }
        }

        std::cout << std::setw(6) << n << " x " << std::setw(4) << n
                  << "  |  " << std::setw(11) << blas_time
                  << "  |  " << std::setw(11) << naive_time
                  << "  |  " << std::setw(5) << std::setprecision(1) << speedup << "x"
                  << "  |  " << (match ? "YES" : "NO") << "\n";
        std::cout << std::setprecision(2);  // 恢复精度
    }
    std::cout << "\n";
}

// ==========================================================================
// NTT 演示
// ==========================================================================

void demo_ntt_kyber()
{
    print_separator("NTT 演示: Kyber 参数 (q=3329, n=256)");

    using Params = ntt::KyberParams;

    // 验证参数正确性
    std::cout << "--- 参数验证 ---\n";
    bool valid = ntt::verify_params<Params>();
    std::cout << "Kyber 参数 (q=" << Params::q << ", n=" << Params::n
              << ", omega=" << Params::omega << "): "
              << (valid ? "✓ 正确" : "✗ 错误") << "\n\n";

    if (!valid) return;

    // 创建 NTT 上下文
    ntt::NTTContext<Params> ctx;

    // --- 1. 基本 NTT 正逆变换验证 ---
    std::cout << "--- 1. NTT 正逆变换验证 ---\n";
    std::cout << "验证: INTT(NTT(a)) == a\n\n";

    // 创建测试多项式（仅显示前 16 个系数）
    std::vector<uint32_t> a(Params::n);
    for (uint32_t i = 0; i < Params::n; ++i) {
        a[i] = i % Params::q;  // 简单的测试数据：0, 1, 2, ..., 255
    }

    // 保存原始数据
    std::vector<uint32_t> a_orig(a);

    // 打印原始多项式（前 16 项）
    std::cout << "原始多项式 a（前 16 项）: ";
    for (int i = 0; i < 16; ++i) std::cout << a[i] << " ";
    std::cout << "...\n";

    // 正向 NTT
    Timer timer;
    timer.start();
    ctx.forward(a.data());
    double fwd_time = timer.elapsed_us();
    std::cout << "NTT 域（前 16 项）: ";
    for (int i = 0; i < 16; ++i) std::cout << a[i] << " ";
    std::cout << "...\n";
    std::cout << "正向 NTT 耗时: " << std::fixed << std::setprecision(1)
              << fwd_time << " μs\n";

    // 逆向 NTT
    timer.start();
    ctx.inverse(a.data());
    double inv_time = timer.elapsed_us();
    std::cout << "恢复后（前 16 项）: ";
    for (int i = 0; i < 16; ++i) std::cout << a[i] << " ";
    std::cout << "...\n";
    std::cout << "逆向 NTT 耗时: " << inv_time << " μs\n";

    // 验证恢复正确性
    bool correct = (a == a_orig);
    std::cout << "INTT(NTT(a)) == a: " << (correct ? "✓ 正确" : "✗ 错误") << "\n\n";

    // --- 2. 多项式乘法 ---
    std::cout << "--- 2. NTT 多项式乘法 ---\n";
    std::cout << "在 Z_q[x]/(x^n + 1) 上计算 c = a * b\n\n";

    // 使用较小的多项式便于对比验证
    const size_t test_n = 16;  // 使用小规模进行可视化演示
    std::vector<uint32_t> pa(test_n), pb(test_n), pc_ntt(test_n), pc_naive(test_n);

    // 初始化测试数据
    for (size_t i = 0; i < test_n; ++i) {
        pa[i] = (i + 1) % Params::q;
        pb[i] = (2 * i + 1) % Params::q;
    }

    std::cout << "a = [";
    for (size_t i = 0; i < test_n; ++i) {
        std::cout << pa[i];
        if (i < test_n - 1) std::cout << ", ";
    }
    std::cout << "]\n";

    std::cout << "b = [";
    for (size_t i = 0; i < test_n; ++i) {
        std::cout << pb[i];
        if (i < test_n - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";

    // NTT 乘法（使用完整 256 维上下文，但只取前 test_n 项用于演示）
    // 这里使用朴素乘法来验证
    ntt::naive_poly_mul(pa.data(), pb.data(), pc_naive.data(), test_n, Params::q);

    std::cout << "c = a * b mod (x^" << test_n << " + 1, " << Params::q << "):\n  [";
    for (size_t i = 0; i < test_n; ++i) {
        std::cout << pc_naive[i];
        if (i < test_n - 1) std::cout << ", ";
    }
    std::cout << "]\n\n";

    // --- 3. 完整多项式乘法性能测试 ---
    std::cout << "--- 3. 多项式乘法性能 (n=" << Params::n << ") ---\n";

    std::vector<uint32_t> poly_a(Params::n), poly_b(Params::n), poly_c(Params::n);
    for (uint32_t i = 0; i < Params::n; ++i) {
        poly_a[i] = (i * 7 + 3) % Params::q;
        poly_b[i] = (i * 13 + 5) % Params::q;
    }

    // NTT 乘法
    const int repeat = 1000;
    timer.start();
    for (int i = 0; i < repeat; ++i) {
        ctx.poly_mul(poly_a.data(), poly_b.data(), poly_c.data());
    }
    double ntt_time = timer.elapsed_ms() / repeat;

    // 朴素乘法（仅运行一次，因为太慢）
    std::vector<uint32_t> poly_c_naive(Params::n);
    timer.start();
    ntt::naive_poly_mul(poly_a.data(), poly_b.data(), poly_c_naive.data(),
                        Params::n, Params::q);
    double naive_time = timer.elapsed_ms();

    std::cout << "NTT 多项式乘法: " << std::fixed << std::setprecision(3)
              << ntt_time << " ms (平均 " << repeat << " 次)\n";
    std::cout << "朴素多项式乘法: " << naive_time << " ms (单次)\n";

    if (naive_time > 0.001) {
        std::cout << "加速比: " << std::setprecision(1) << naive_time / ntt_time << "x\n";
    }

    // 验证结果一致性
    bool mul_correct = (poly_c == poly_c_naive);
    std::cout << "结果一致性: " << (mul_correct ? "✓ 正确" : "✗ 错误") << "\n\n";
}

void demo_ntt_dilithium()
{
    print_separator("NTT 演示: Dilithium 参数 (q=8380417, n=256)");

    using Params = ntt::DilithiumParams;

    // 验证参数
    std::cout << "--- 参数验证 ---\n";
    bool valid = ntt::verify_params<Params>();
    std::cout << "Dilithium 参数 (q=" << Params::q << ", n=" << Params::n
              << ", omega=" << Params::omega << "): "
              << (valid ? "✓ 正确" : "✗ 错误") << "\n\n";

    if (!valid) return;

    // 创建 NTT 上下文
    ntt::NTTContext<Params> ctx;

    // 简单的正逆变换验证
    std::cout << "--- NTT 正逆变换验证 ---\n";
    std::vector<uint32_t> a(Params::n);
    for (uint32_t i = 0; i < Params::n; ++i) {
        a[i] = (i * 31 + 17) % Params::q;
    }
    std::vector<uint32_t> a_orig(a);

    ctx.forward(a.data());
    ctx.inverse(a.data());

    bool correct = (a == a_orig);
    std::cout << "INTT(NTT(a)) == a: " << (correct ? "✓ 正确" : "✗ 错误") << "\n\n";

    // 性能测试
    std::cout << "--- 多项式乘法性能 ---\n";
    std::vector<uint32_t> poly_a(Params::n), poly_b(Params::n), poly_c(Params::n);
    for (uint32_t i = 0; i < Params::n; ++i) {
        poly_a[i] = (i * 7 + 3) % Params::q;
        poly_b[i] = (i * 13 + 5) % Params::q;
    }

    const int repeat = 1000;
    Timer timer;
    timer.start();
    for (int i = 0; i < repeat; ++i) {
        ctx.poly_mul(poly_a.data(), poly_b.data(), poly_c.data());
    }
    double ntt_time = timer.elapsed_ms() / repeat;

    std::cout << "NTT 多项式乘法: " << std::fixed << std::setprecision(3)
              << ntt_time << " ms (平均 " << repeat << " 次)\n\n";
}

// ==========================================================================
// OpenBLAS 环境信息
// ==========================================================================

void print_openblas_info()
{
    print_separator("OpenBLAS 环境信息");

    // 获取 OpenBLAS 版本信息
    std::cout << "OpenBLAS 版本: " << openblas_get_config() << "\n";
    std::cout << "CPU 内核数: " << openblas_get_num_procs() << "\n";
    std::cout << "当前线程数: " << openblas_get_num_threads() << "\n";
    std::cout << "CPU 型号: " << openblas_get_corename() << "\n\n";
}

// ==========================================================================
// 主函数
// ==========================================================================

int main(int argc, char* argv[])
{
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║            OpenBLAS 演示项目 — BLAS 基础算子与 NTT            ║\n";
    std::cout << "║           面向密码学研究的线性代数与多项式运算演示             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";

    // 解析命令行参数
    std::string mode = "all";
    if (argc > 1) {
        mode = argv[1];
    }

    // 显示 OpenBLAS 环境信息
    print_openblas_info();

    // 根据参数运行不同的演示模块
    if (mode == "all" || mode == "--blas") {
        demo_blas_level1();
        demo_blas_level2();
        demo_blas_level3();
    }

    if (mode == "all" || mode == "--bench") {
        benchmark_gemm();
    }

    if (mode == "all" || mode == "--ntt") {
        demo_ntt_kyber();
        demo_ntt_dilithium();
    }

    std::cout << "\n演示完成！\n";
    return 0;
}
