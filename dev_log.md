# Init (2026/6/17 16:00)
  ✦ 项目已创建完成，结构如下：

      1 OpenBLAS/
      2 ├── CMakeLists.txt        # CMake 构建系统（自动查找 OpenBLAS）
      3 ├── README.md             # 完整文档（编译、运行、调试、FAQ）
      4 ├── .gitignore
      5 ├── include/
      6 │   ├── blas_demo.h       # BLAS Level 1/2/3 接口（8 个函数）
      7 │   └── ntt.h             # NTT 接口（Kyber + Dilithium 参数）
      8 └── src/
      9     ├── main.cpp           # 主程序（支持 --blas / --ntt / --bench 参数）
     10     ├── blas_demo.cpp      # CBLAS 封装，每个调用附详细中文注释
     11     └── ntt.cpp            # 迭代 CT/GS 蝶形 NTT 实现

    核心内容：
     - BLAS: DDOT、DAXPY、DSCAL、DNRM2、DGEMV、DGEMM、DSYMM，均通过 cblas_* 接口调用 OpenBLAS
     - NTT: 支持 Kyber (q=3329, n=256, ω=17) 和 Dilithium (q=8380417) 两组参数，包含 Cooley-Tukey 正向变换、Gentleman-Sande
       逆变换、逐点乘法、完整多项式乘法流程，以及朴素 O(n²) 验证实现
     - 调试: README 中包含 GDB、Valgrind、AddressSanitizer、perf 的详细使用说明

    编译命令：

     1 mkdir -p build && cd build
     2 cmake .. -DCMAKE_BUILD_TYPE=Release
     3 make -j$(nproc)
     4 ./openblas_demo

---

## Fix (2026/6/17) — 修复 `CBLAS_INT` 类型未定义的编译错误

### 问题描述

编译 `src/blas_demo.cpp` 时报错：

```
error: 'CBLAS_INT' does not name a type; did you mean 'CBLAS_SIDE'?
```

所有调用 CBLAS 函数的地方（`cblas_ddot`、`cblas_daxpy`、`cblas_dscal`、`cblas_dnrm2`、
`cblas_dasum`、`cblas_dgemv`、`cblas_dgemm`、`cblas_dsymm`）中，`static_cast<CBLAS_INT>(...)` 均报错，
共涉及 **19 处**类型转换。

### 原因分析

`CBLAS_INT` 是 OpenBLAS 较新版本（≥ 0.3.17）引入的类型别名，用于统一 CBLAS 接口的整型参数。
但不同版本的 OpenBLAS 对 CBLAS 整型参数的定义不一致：

| OpenBLAS 版本   | 整型类型    | 说明                                                       |
|----------------|------------|-----------------------------------------------------------|
| < 0.3.0        | `int`      | 原始 BLAS 接口，使用 C 标准 `int`                            |
| 0.3.0 ~ 0.3.16 | `blasint`  | OpenBLAS 自定义类型（通常为 `int`，INTERFACE64 模式下为 `int64_t`） |
| ≥ 0.3.17       | `CBLAS_INT`| 新的标准别名，等价于 `blasint`                                 |

用户系统的 OpenBLAS 版本低于 0.3.17，头文件中未定义 `CBLAS_INT`，导致编译失败。

### 修复方案

将所有 `static_cast<CBLAS_INT>(...)` 替换为 `static_cast<int>(...)`。

**选择 `int` 的理由：**
1. CBLAS 标准（Netlib reference）定义的参数类型就是 `int`
2. `int` 在所有版本的 OpenBLAS 中都是合法的参数类型
   （即使内部 typedef 为 `blasint`，`blasint` 默认就是 `int`）
3. 兼容性最好，覆盖所有 OpenBLAS 版本

**潜在风险：**
- 当矩阵维度超过 `INT_MAX`（约 21 亿）时，`int` 会溢出。
  但对于本演示项目的规模（最大 2048×2048）完全不存在此问题
- 如需支持超大矩阵（ILP64 模式），应在 CMakeLists.txt 中检测 INTERFACE64 宏并相应调整类型

### 修改文件

| 文件                | 修改内容                                        |
|--------------------|-------------------------------------------------|
| `src/blas_demo.cpp` | 全局替换 `static_cast<CBLAS_INT>` → `static_cast<int>`（19 处） |

### 验证

修复后重新编译：

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

---

## Fix (2026/6/17) — 修复 NTT 实现：从 Cyclic NTT 改为 Negacyclic NTT

### 问题描述

运行程序后发现两个错误：

1. **Kyber NTT 多项式乘法结果不一致**：
   ```
   NTT 多项式乘法: 0.033 ms (平均 1000 次)
   朴素多项式乘法: 0.449 ms (单次)
   结果一致性: ✗ 错误
   ```
   `poly_mul`（NTT 加速）和 `naive_poly_mul`（朴素 O(n²)）的结果不匹配。

2. **Dilithium 参数验证失败**：
   ```
   错误: omega^n = 8380416 != 1 (mod 8380417)
   Dilithium 参数 (q=8380417, n=256, omega=1753): ✗ 错误
   ```
   `verify_params` 要求 `omega^n ≡ 1`，但 Dilithium 的 `1753^256 ≡ -1 (mod 8380417)`。

### 原因分析（核心密码学问题）

原始实现使用的是 **标准（cyclic）NTT**，它计算的是环 Z_q[x]/(x^n - 1) 上的
**循环卷积**。但格密码方案（Kyber、Dilithium）工作在环 Z_q[x]/(x^n + 1) 上，
需要的是 **negacyclic 卷积**（反循环卷积）。

两者的关键区别：

| 性质              | Cyclic NTT             | Negacyclic NTT              |
|-------------------|------------------------|------------------------------|
| 约减多项式         | x^n - 1               | x^n + 1                      |
| 所需单位根         | ω: ω^n = 1            | ψ: ψ^n = -1（2n 次单位根）   |
| 蝶形旋转因子       | ω^k                   | ψ 的奇数次幂                 |
| 逐点相乘对应的卷积  | cyclic (mod x^n - 1)  | negacyclic (mod x^n + 1)     |

原始代码使用 ω^n = 1 的标准 NTT，其逐点相乘等价于 mod (x^n - 1) 的循环卷积，
而非 `naive_poly_mul` 中 mod (x^n + 1) 的 negacyclic 卷积。因此两者结果不一致。

对于 Dilithium，`omega = 1753` 满足 `omega^256 = -1`，它是 2n 次（512 次）
单位根而非 n 次单位根。原始的 `verify_params` 检查 `omega^n == 1` 自然失败。

### 修复方案：Pre/Post-Multiplication 方法

采用经典的 **negacyclic NTT 实现方法**：

设 ψ 为 2n 次本原单位根（ψ^n ≡ -1 mod q），ω = ψ^2 为 n 次单位根。

**正向 Negacyclic NTT**：
1. 预乘：`poly[i] *= ψ^i mod q`（negacyclic twist）
2. 标准 NTT（使用 ω = ψ^2 作为旋转因子）

**逆向 Negacyclic NTT**：
1. 标准 INTT（使用 ω^{-1}）
2. 后乘：`poly[i] *= ψ^{-i} mod q`（撤销 twist）

**多项式乘法**：
```
c = INTT_negacyclic(NTT_negacyclic(a) ⊙ NTT_negacyclic(b))
  = a * b in Z_q[x]/(x^n + 1)
```

**ψ 的推导**（两种情况）：
- **Kyber**（omega=17, omega^256 = 1）：psi = sqrt(17) mod 3329，使用 Tonelli-Shanks 算法
- **Dilithium**（omega=1753, omega^256 = -1）：psi = omega = 1753（直接可用）

### 修改文件

| 文件            | 修改内容                                                            |
|----------------|---------------------------------------------------------------------|
| `include/ntt.h` | 添加 `mod_sqrt` 声明；更新 NTTContext 成员变量（ψ, ω_ntt, ψ_inv）；更新文档注释 |
| `src/ntt.cpp`   | 添加 `mod_sqrt`（Tonelli-Shanks）实现；重写构造函数推导 ψ；forward 添加预乘；inverse 添加后乘；重写 verify_params |

### 具体代码变更

1. **新增 `mod_sqrt` 函数**：Tonelli-Shanks 算法，O(log² q) 复杂度
2. **构造函数**：自动检测 omega 是 n 次还是 2n 次单位根，推导 ψ 和 ω_ntt
3. **`forward()`**：在比特反转和蝶形运算前，增加 `poly[i] *= ψ^i` 预乘步骤
4. **`inverse()`**：在蝶形运算和归一化后，增加 `poly[i] *= ψ^{-i}` 后乘步骤
5. **`verify_params()`**：兼容两种 omega 类型，验证 ψ^n ≡ -1 (mod q)

### 验证

修复后重新编译运行：

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)
./openblas_demo --ntt
```

预期结果：
- Kyber 参数验证：✓ 正确
- INTT(NTT(a)) == a：✓ 正确
- 结果一致性（NTT vs 朴素）：✓ 正确
- Dilithium 参数验证：✓ 正确
