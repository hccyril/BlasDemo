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
