/**
 * @file blas_demo.h
 * @brief BLAS 基础算子演示接口
 *
 * 本文件封装了 OpenBLAS (CBLAS) 的核心函数调用，覆盖 BLAS Level 1/2/3：
 *   - Level 1: 向量-向量操作（点积、缩放、范数等）
 *   - Level 2: 矩阵-向量操作（矩阵向量乘法）
 *   - Level 3: 矩阵-矩阵操作（通用矩阵乘法 DGEMM）
 *
 * 对于密码学研究：
 *   - 格密码（Lattice-based Cryptography）中大量使用矩阵运算
 *   - 同态加密（FHE）方案的 bootstrapping 涉及大矩阵乘法
 *   - BLAS Level 3 (GEMM) 是性能关键路径
 *
 * @note 所有矩阵采用列优先（Column-major）存储，与 Fortran/LAPACK 一致
 *       如需行优先，请在调用时设置 CblasRowMajor
 */

#ifndef BLAS_DEMO_H
#define BLAS_DEMO_H

#include <cstddef>
#include <vector>
#include <string>

namespace blas_demo {

// ==========================================================================
// 工具函数
// ==========================================================================

/**
 * @brief 打印矩阵内容（用于调试和演示）
 * @param name   矩阵名称（显示在输出前）
 * @param data   矩阵数据指针（列优先存储）
 * @param rows   行数
 * @param cols   列数
 * @param lda    主维度（leading dimension），列优先时 lda >= rows
 */
void print_matrix(const std::string& name, const double* data,
                  size_t rows, size_t cols, size_t lda);

/**
 * @brief 打印向量内容
 * @param name   向量名称
 * @param data   向量数据指针
 * @param n      向量长度
 */
void print_vector(const std::string& name, const double* data, size_t n);

/**
 * @brief 随机初始化向量（均匀分布在 [lo, hi]）
 * @param data  输出向量
 * @param n     向量长度
 * @param lo    下界（默认 0.0）
 * @param hi    上界（默认 1.0）
 */
void random_vector(double* data, size_t n, double lo = 0.0, double hi = 1.0);

/**
 * @brief 随机初始化矩阵（均匀分布在 [lo, hi]）
 * @param data  输出矩阵（列优先存储）
 * @param rows  行数
 * @param cols  列数
 * @param lo    下界（默认 -1.0）
 * @param hi    上界（默认 1.0）
 */
void random_matrix(double* data, size_t rows, size_t cols,
                   double lo = -1.0, double hi = 1.0);

// ==========================================================================
// BLAS Level 1: 向量-向量操作
// ==========================================================================

/**
 * @brief 向量点积（内积）: result = x^T * y
 *
 * 对应 BLAS 函数: ddot
 * 计算复杂度: O(n)
 *
 * 在密码学中的应用：
 *   - 格上的内积计算（Inner Product Argument）
 *   - LWE 问题中的 <s, a> 计算
 *
 * @param x     第一个向量
 * @param y     第二个向量
 * @param n     向量长度
 * @return      点积结果
 */
double dot_product(const double* x, const double* y, size_t n);

/**
 * @brief 向量缩放: y = alpha * x + y (AXPY 操作)
 *
 * 对应 BLAS 函数: daxpy
 * 计算复杂度: O(n)
 *
 * 这是最基本的 BLAS 操作之一，广泛用于：
 *   - 梯度下降等优化算法
 *   - 线性组合计算
 *
 * @param alpha  缩放因子
 * @param x      输入向量（不被修改）
 * @param y      输入/输出向量（原地更新）
 * @param n      向量长度
 */
void axpy(double alpha, const double* x, double* y, size_t n);

/**
 * @brief 向量缩放: x = alpha * x (SCAL 操作)
 *
 * 对应 BLAS 函数: dscal
 * 计算复杂度: O(n)
 *
 * @param alpha  缩放因子
 * @param x      输入/输出向量（原地缩放）
 * @param n      向量长度
 */
void scal(double alpha, double* x, size_t n);

/**
 * @brief 向量欧几里得范数（L2 范数）: ||x||_2
 *
 * 对应 BLAS 函数: dnrm2
 * 计算复杂度: O(n)
 *
 * 在密码学中的应用：
 *   - 格基约化（LLL, BKZ）中向量的长度计算
 *   - 高斯采样中的范数约束
 *
 * @param x  输入向量
 * @param n  向量长度
 * @return   L2 范数值
 */
double vector_norm(const double* x, size_t n);

/**
 * @brief 向量元素绝对值之和（L1 范数）: sum(|x_i|)
 *
 * 对应 BLAS 函数: dasum
 * 计算复杂度: O(n)
 *
 * @param x  输入向量
 * @param n  向量长度
 * @return   L1 范数值
 */
double vector_asum(const double* x, size_t n);

// ==========================================================================
// BLAS Level 2: 矩阵-向量操作
// ==========================================================================

/**
 * @brief 通用矩阵-向量乘法: y = alpha * A * x + beta * y (GEMV)
 *
 * 对应 BLAS 函数: dgemv
 * 计算复杂度: O(m * n)
 *
 * 在密码学中的应用：
 *   - LWE 加密: b = A * s + e (mod q)
 *   - 签名方案中的矩阵向量乘法
 *   - 零知识证明中的线性约束验证
 *
 * @param A      矩阵（列优先，m x n）
 * @param x      输入向量（长度 n）
 * @param y      输出向量（长度 m）
 * @param m      矩阵行数
 * @param n      矩阵列数
 * @param lda    矩阵 A 的主维度
 * @param alpha  缩放因子 alpha（默认 1.0）
 * @param beta   缩放因子 beta（默认 0.0）
 */
void gemv(const double* A, const double* x, double* y,
          size_t m, size_t n, size_t lda,
          double alpha = 1.0, double beta = 0.0);

// ==========================================================================
// BLAS Level 3: 矩阵-矩阵操作
// ==========================================================================

/**
 * @brief 通用矩阵乘法: C = alpha * A * B + beta * C (GEMM)
 *
 * 对应 BLAS 函数: dgemm
 * 计算复杂度: O(m * n * k)
 *
 * 这是 BLAS 中最重要的函数，也是 OpenBLAS 优化的核心。
 * OpenBLAS 通过以下技术优化 GEMM：
 *   - 分块（Blocking/Tiling）适应 CPU 缓存层次
 *   - SIMD 向量化（AVX2, AVX-512, NEON 等）
 *   - 多线程并行
 *   - 汇编级微内核（Micro-kernel）
 *
 * 在密码学中的应用：
 *   - 多线性映射（Multilinear Maps）
 *   - 基于格的签名方案批量验证
 *   - FHE 方案中的矩阵乘法（如 BGV, CKKS）
 *
 * @param A      左矩阵（m x k，列优先）
 * @param B      右矩阵（k x n，列优先）
 * @param C      输出矩阵（m x n，列优先）
 * @param m      A 的行数 / C 的行数
 * @param n      B 的列数 / C 的列数
 * @param k      A 的列数 / B 的行数
 * @param lda    A 的主维度
 * @param ldb    B 的主维度
 * @param ldc    C 的主维度
 * @param alpha  缩放因子 alpha（默认 1.0）
 * @param beta   缩放因子 beta（默认 0.0）
 */
void gemm(const double* A, const double* B, double* C,
          size_t m, size_t n, size_t k,
          size_t lda, size_t ldb, size_t ldc,
          double alpha = 1.0, double beta = 0.0);

/**
 * @brief 对称矩阵乘法: C = alpha * A * B + beta * C，其中 A 为对称矩阵
 *
 * 对应 BLAS 函数: dsymm
 * 计算复杂度: O(m * n * k)，但利用对称性减少约一半的内存访问
 *
 * 在密码学中，Gram 矩阵（G = B^T * B）是对称的，
 * 格基约化算法中频繁使用对称矩阵乘法。
 *
 * @param A      对称矩阵（m x m，列优先，仅使用上三角或下三角）
 * @param B      一般矩阵（m x n，列优先）
 * @param C      输出矩阵（m x n，列优先）
 * @param m      A 的维度 / C 的行数
 * @param n      B 的列数 / C 的列数
 * @param lda    A 的主维度
 * @param ldb    B 的主维度
 * @param ldc    C 的主维度
 * @param upper  true 表示使用 A 的上三角部分，false 表示下三角
 * @param alpha  缩放因子（默认 1.0）
 * @param beta   缩放因子（默认 0.0）
 */
void symm(const double* A, const double* B, double* C,
          size_t m, size_t n,
          size_t lda, size_t ldb, size_t ldc,
          bool upper = true,
          double alpha = 1.0, double beta = 0.0);

// ==========================================================================
// 朴素实现（不使用 BLAS，用于性能对比）
// ==========================================================================

/**
 * @brief 朴素矩阵乘法: C = A * B（不使用 BLAS，三层循环实现）
 *
 * 用于与 cblas_dgemm 进行性能对比，展示 OpenBLAS 优化带来的加速效果。
 * 矩阵采用列优先存储，与 BLAS 接口一致。
 *
 * @param A    左矩阵（m x k，列优先）
 * @param B    右矩阵（k x n，列优先）
 * @param C    输出矩阵（m x n，列优先）
 * @param m    A 的行数 / C 的行数
 * @param n    B 的列数 / C 的列数
 * @param k    A 的列数 / B 的行数
 * @param lda  A 的主维度
 * @param ldb  B 的主维度
 * @param ldc  C 的主维度
 */
void naive_gemm(const double* A, const double* B, double* C,
                size_t m, size_t n, size_t k,
                size_t lda, size_t ldb, size_t ldc);

} // namespace blas_demo

#endif // BLAS_DEMO_H
