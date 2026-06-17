/**
 * @file blas_demo.cpp
 * @brief BLAS 基础算子实现 — 封装 OpenBLAS (CBLAS) 调用
 *
 * 本文件演示如何通过 CBLAS C 接口调用 OpenBLAS 的高性能线性代数例程。
 * 每个函数对应一个标准的 BLAS 操作，并添加了详细的中文注释说明
 * 参数含义、性能特征以及在密码学中的应用场景。
 */

#include "blas_demo.h"

// CBLAS 头文件 —— OpenBLAS 提供的 C 语言接口
// cblas.h 定义了所有 BLAS 函数的 C 接口，命名规则为 cblas_<精度><操作>
// 例如: cblas_ddot = double precision dot product
extern "C" {
#include <cblas.h>
}

#include <iostream>
#include <iomanip>
#include <random>
#include <cstring>

namespace blas_demo {

// ==========================================================================
// 工具函数实现
// ==========================================================================

void print_matrix(const std::string& name, const double* data,
                  size_t rows, size_t cols, size_t lda)
{
    std::cout << name << " (" << rows << " x " << cols << "):\n";
    std::cout << std::fixed << std::setprecision(4);

    for (size_t i = 0; i < rows; ++i) {
        std::cout << "  [";
        for (size_t j = 0; j < cols; ++j) {
            // 列优先存储：元素 (i, j) 位于 data[j * lda + i]
            std::cout << std::setw(10) << data[j * lda + i];
            if (j < cols - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << std::endl;
}

void print_vector(const std::string& name, const double* data, size_t n)
{
    std::cout << name << " (长度 " << n << "): [";
    std::cout << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < n; ++i) {
        std::cout << data[i];
        if (i < n - 1) std::cout << ", ";
    }
    std::cout << "]\n" << std::endl;
}

void random_vector(double* data, size_t n, double lo, double hi)
{
    // 使用 Mersenne Twister 生成高质量随机数
    static std::mt19937_64 rng(42);  // 固定种子以便复现
    std::uniform_real_distribution<double> dist(lo, hi);
    for (size_t i = 0; i < n; ++i) {
        data[i] = dist(rng);
    }
}

void random_matrix(double* data, size_t rows, size_t cols, double lo, double hi)
{
    static std::mt19937_64 rng(123);
    std::uniform_real_distribution<double> dist(lo, hi);
    // 列优先存储：按列填充
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            data[j * rows + i] = dist(rng);
        }
    }
}

// ==========================================================================
// BLAS Level 1: 向量-向量操作
// ==========================================================================

double dot_product(const double* x, const double* y, size_t n)
{
    /**
     * cblas_ddot 参数说明：
     *   N    - 向量长度
     *   X    - 第一个向量指针
     *   incX - x 的步长（stride），1 表示连续存储
     *   Y    - 第二个向量指针
     *   incY - y 的步长
     *
     * OpenBLAS 内部优化：
     *   - 对于大步长（incX > 1）有专门的 kernel
     *   - 使用 SIMD 指令（如 AVX2 的 vmulpd + vaddpd）进行向量化
     */
    return cblas_ddot(
        static_cast<CBLAS_INT>(n),  // N
        x, 1,                       // X, incX
        y, 1                        // Y, incY
    );
}

void axpy(double alpha, const double* x, double* y, size_t n)
{
    /**
     * cblas_daxpy: y := alpha * x + y
     *
     * 参数说明：
     *   N     - 向量长度
     *   alpha - 缩放因子
     *   X     - 输入向量
     *   incX  - x 的步长
     *   Y     - 输入/输出向量（原地修改）
     *   incY  - y 的步长
     *
     * 这是 BLAS 中使用最频繁的操作之一。
     * OpenBLAS 的 daxpy kernel 通常可以达到内存带宽的理论峰值。
     */
    cblas_daxpy(
        static_cast<CBLAS_INT>(n),  // N
        alpha,                      // alpha
        x, 1,                       // X, incX
        y, 1                        // Y, incY
    );
}

void scal(double alpha, double* x, size_t n)
{
    /**
     * cblas_dscal: x := alpha * x
     *
     * 参数说明：
     *   N     - 向量长度
     *   alpha - 缩放因子
     *   X     - 输入/输出向量（原地修改）
     *   incX  - x 的步长
     */
    cblas_dscal(
        static_cast<CBLAS_INT>(n),  // N
        alpha,                      // alpha
        x, 1                        // X, incX
    );
}

double vector_norm(const double* x, size_t n)
{
    /**
     * cblas_dnrm2: ||x||_2 = sqrt(sum(x_i^2))
     *
     * 参数说明同 ddot。
     *
     * 实现细节：OpenBLAS 使用数值稳定的算法避免溢出/下溢，
     * 类似于 LAPACK 的 dlapy2 中的缩放策略。
     */
    return cblas_dnrm2(
        static_cast<CBLAS_INT>(n),  // N
        x, 1                        // X, incX
    );
}

double vector_asum(const double* x, size_t n)
{
    /**
     * cblas_dasum: sum(|x_i|)
     *
     * L1 范数计算，在稀疏优化（LASSO）和压缩感知中有广泛应用。
     */
    return cblas_dasum(
        static_cast<CBLAS_INT>(n),  // N
        x, 1                        // X, incX
    );
}

// ==========================================================================
// BLAS Level 2: 矩阵-向量操作
// ==========================================================================

void gemv(const double* A, const double* x, double* y,
          size_t m, size_t n, size_t lda,
          double alpha, double beta)
{
    /**
     * cblas_dgemv: y := alpha * op(A) * x + beta * y
     *
     * 参数说明：
     *   Order  - 存储顺序（CblasColMajor = 列优先 = Fortran 顺序）
     *   TransA - 是否转置 A：
     *              CblasNoTrans   -> y = alpha * A * x + beta * y
     *              CblasTrans     -> y = alpha * A^T * x + beta * y
     *   M      - A 的行数
     *   N      - A 的列数
     *   alpha  - 缩放因子
     *   A      - 矩阵指针
     *   lda    - A 的主维度（列优先时 lda >= M）
     *   X      - 输入向量
     *   incX   - x 的步长
     *   beta   - y 的缩放因子（beta=0 时 y 的初值被忽略）
     *   Y      - 输出向量
     *   incY   - y 的步长
     *
     * 在密码学 LWE 问题中，核心运算 b = As + e 就是 GEMV 操作。
     * OpenBLAS 对 GEMV 的优化包括：
     *   - 行/列访问模式自动检测与优化
     *   - 多线程并行（当向量足够大时）
     */
    cblas_dgemv(
        CblasColMajor,                                  // 列优先存储
        CblasNoTrans,                                   // 不转置
        static_cast<CBLAS_INT>(m),                      // M
        static_cast<CBLAS_INT>(n),                      // N
        alpha,                                          // alpha
        A, static_cast<CBLAS_INT>(lda),                 // A, lda
        x, 1,                                           // X, incX
        beta,                                           // beta
        y, 1                                            // Y, incY
    );
}

// ==========================================================================
// BLAS Level 3: 矩阵-矩阵操作
// ==========================================================================

void gemm(const double* A, const double* B, double* C,
          size_t m, size_t n, size_t k,
          size_t lda, size_t ldb, size_t ldc,
          double alpha, double beta)
{
    /**
     * cblas_dgemm: C := alpha * op(A) * op(B) + beta * C
     *
     * 这是 BLAS 中最核心、优化最充分的函数。
     * OpenBLAS 的 DGEMM 性能通常可以达到 CPU 理论峰值的 80-95%。
     *
     * 参数说明：
     *   Order  - 存储顺序
     *   TransA - A 的转置选项
     *   TransB - B 的转置选项
     *   M      - op(A) 的行数 / C 的行数
     *   N      - op(B) 的列数 / C 的列数
     *   K      - op(A) 的列数 / op(B) 的行数
     *   alpha  - 缩放因子
     *   A      - 左矩阵
     *   lda    - A 的主维度
     *   B      - 右矩阵
     *   ldb    - B 的主维度
     *   beta   - C 的缩放因子
     *   C      - 输出矩阵
     *   ldc    - C 的主维度
     *
     * OpenBLAS DGEMM 内部实现（简化）：
     *   1. 将矩阵分块（Blocking），适应 L1/L2/L3 缓存
     *   2. 最内层使用汇编编写的微内核（Micro-kernel）
     *      - 例如 Haswell 上使用 AVX2 + FMA，每次处理 8x8 子块
     *   3. 多线程：将大矩阵切分为子任务分配给各线程
     *   4. 内存打包（Packing）：将分块数据按访问顺序排列到连续内存
     *
     * 性能参考（Intel Xeon, AVX-512）：
     *   - n=1000: ~100 GFLOPS（单核）, ~800 GFLOPS（多核）
     *   - n=4000: 接近理论峰值
     */
    cblas_dgemm(
        CblasColMajor,                                  // 列优先
        CblasNoTrans,                                   // A 不转置
        CblasNoTrans,                                   // B 不转置
        static_cast<CBLAS_INT>(m),                      // M
        static_cast<CBLAS_INT>(n),                      // N
        static_cast<CBLAS_INT>(k),                      // K
        alpha,                                          // alpha
        A, static_cast<CBLAS_INT>(lda),                 // A, lda
        B, static_cast<CBLAS_INT>(ldb),                 // B, ldb
        beta,                                           // beta
        C, static_cast<CBLAS_INT>(ldc)                  // C, ldc
    );
}

void symm(const double* A, const double* B, double* C,
          size_t m, size_t n,
          size_t lda, size_t ldb, size_t ldc,
          bool upper,
          double alpha, double beta)
{
    /**
     * cblas_dsymm: C := alpha * A * B + beta * C，其中 A 为对称矩阵
     *
     * 参数说明：
     *   Side   - 对称矩阵在哪一侧：
     *              CblasLeft  -> C = alpha * A * B + beta * C（A 在左）
     *              CblasRight -> C = alpha * B * A + beta * C（A 在右）
     *   Uplo   - 使用上三角还是下三角：
     *              CblasUpper -> 仅访问 A 的上三角部分
     *              CblasLower -> 仅访问 A 的下三角部分
     *
     * 对称矩阵乘法在格密码中的应用：
     *   - Gram 矩阵 G = B^T * B 是对称正定的
     *   - LLL/BKZ 算法中频繁使用 G * v 的计算
     *   - 利用对称性可减少约一半的内存访问和计算
     */
    cblas_dsymm(
        CblasColMajor,                                  // 列优先
        CblasLeft,                                      // A 在左侧
        upper ? CblasUpper : CblasLower,                // 上/下三角
        static_cast<CBLAS_INT>(m),                      // M
        static_cast<CBLAS_INT>(n),                      // N
        alpha,                                          // alpha
        A, static_cast<CBLAS_INT>(lda),                 // A, lda
        B, static_cast<CBLAS_INT>(ldb),                 // B, ldb
        beta,                                           // beta
        C, static_cast<CBLAS_INT>(ldc)                  // C, ldc
    );
}

} // namespace blas_demo
