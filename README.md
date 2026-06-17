# OpenBLAS 演示项目

面向密码学研究的 **BLAS 基础算子** 与 **数论变换（NTT）** 演示项目。

## 项目概述

本项目演示了以下核心内容：

| 模块 | 内容 | 密码学应用 |
|------|------|-----------|
| BLAS Level 1 | 向量点积、AXPY、范数 | LWE 内积、格基约化 |
| BLAS Level 2 | 矩阵-向量乘法 (GEMV) | LWE 加密 b=As+e |
| BLAS Level 3 | 矩阵乘法 (GEMM) | FHE 矩阵运算、批量签名验证 |
| NTT | 数论变换 + 多项式乘法 | Kyber/ML-KEM、Dilithium/ML-DSA |

### 技术栈

- **语言**: C++17
- **线性代数库**: OpenBLAS（通过 CBLAS C 接口调用）
- **构建系统**: CMake 3.14+
- **NTT**: 自行实现的教学级 NTT（支持 Kyber 和 Dilithium 参数）

---

## 目录结构

```
OpenBLAS/
├── CMakeLists.txt          # CMake 构建配置
├── README.md               # 本文件
├── include/                # 头文件
│   ├── blas_demo.h         # BLAS 演示接口声明
│   └── ntt.h               # NTT 接口声明
├── src/                    # 源代码
│   ├── main.cpp            # 主程序入口
│   ├── blas_demo.cpp       # BLAS 操作实现（封装 CBLAS 调用）
│   └── ntt.cpp             # NTT 算法实现
├── lib/                    # 自定义库文件（可选，用于放置本地编译的 OpenBLAS）
└── build/                  # 构建输出目录（CMake 生成）
```

---

## 环境依赖

### 1. 编译器

需要支持 C++17 的编译器：

- **GCC** >= 7.0（推荐 >= 9.0）
- **Clang** >= 5.0
- **MSVC** >= 2017（Windows，未测试）

```bash
# 检查 GCC 版本
g++ --version

# 检查 Clang 版本
clang++ --version
```

### 2. OpenBLAS

#### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install libopenblas-dev
```

#### CentOS / RHEL / Fedora

```bash
sudo yum install openblas-devel
# 或 Fedora:
sudo dnf install openblas-devel
```

#### macOS (Homebrew)

```bash
brew install openblas
```

#### 从源码编译（推荐用于性能调优）

```bash
git clone https://github.com/OpenMathLib/OpenBLAS.git
cd OpenBLAS
make -j$(nproc)
sudo make install PREFIX=/opt/OpenBLAS
```

验证安装：

```bash
# 检查头文件是否存在
ls /usr/include/cblas.h
# 或
ls /usr/include/openblas/cblas.h

# 检查库文件
ldconfig -p | grep openblas
```

### 3. CMake

```bash
# Ubuntu/Debian
sudo apt-get install cmake

# macOS
brew install cmake

# 验证版本（需要 >= 3.14）
cmake --version
```

---

## 编译

### 标准编译

```bash
# 进入项目根目录
cd /path/to/OpenBLAS

# 创建构建目录
mkdir -p build && cd build

# 配置 CMake（Release 模式）
cmake .. -DCMAKE_BUILD_TYPE=Release

# 编译
make -j$(nproc)
```

编译成功后，可执行文件位于 `build/openblas_demo`。

### Debug 编译

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

Debug 模式会启用 `-g -O0` 编译选项，便于调试。

### 指定 OpenBLAS 安装路径

如果 OpenBLAS 安装在非标准路径：

```bash
# 方法一：设置环境变量
export OPENBLAS_DIR=/opt/OpenBLAS
cmake .. -DCMAKE_BUILD_TYPE=Release

# 方法二：直接指定库和头文件路径
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DOPENBLAS_LIBRARIES=/opt/OpenBLAS/lib/libopenblas.so \
    -DOPENBLAS_INCLUDE_DIRS=/opt/OpenBLAS/include
```

### 使用 Make 快捷方式

```bash
# 在 build 目录下也可以使用 make 的常用命令
cd build
make -j$(nproc)           # 编译
make clean                 # 清理编译产物
make install               # 安装到系统路径
```

---

## 运行

### 运行全部演示

```bash
./build/openblas_demo
```

### 按模块运行

```bash
# 仅运行 BLAS 演示（Level 1/2/3）
./build/openblas_demo --blas

# 仅运行 NTT 演示（Kyber + Dilithium 参数）
./build/openblas_demo --ntt

# 仅运行矩阵乘法性能基准测试
./build/openblas_demo --bench
```

### 输出示例

程序运行后会输出：

1. **OpenBLAS 环境信息** — 版本号、CPU 内核数、线程数、CPU 型号
2. **BLAS Level 1** — 向量点积、AXPY、SCAL、范数的计算结果
3. **BLAS Level 2** — GEMV 矩阵-向量乘法结果
4. **BLAS Level 3** — GEMM 矩阵乘法、SYMM 对称矩阵乘法结果
5. **性能基准** — 不同规模矩阵乘法的耗时和 GFLOPS
6. **NTT 演示** — 参数验证、正逆变换正确性、多项式乘法性能

---

## 调试

### 使用 GDB 调试

```bash
# 确保使用 Debug 模式编译
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# 启动 GDB
gdb ./openblas_demo

# 常用 GDB 命令：
# (gdb) break main.cpp:230        # 在 main.cpp 第 230 行设置断点
# (gdb) break blas_demo.cpp:85    # 在 BLAS 实现中设置断点
# (gdb) break ntt.cpp:120         # 在 NTT 实现中设置断点
# (gdb) run                       # 运行程序
# (gdb) run --ntt                 # 带参数运行
# (gdb) next                      # 单步执行
# (gdb) step                      # 进入函数
# (gdb) print poly[0]@16          # 打印数组前 16 个元素
# (gdb) info locals               # 查看局部变量
# (gdb) backtrace                 # 查看调用栈
```

### 使用 Valgrind 检查内存问题

```bash
# 安装 Valgrind（如未安装）
sudo apt-get install valgrind

# 检查内存泄漏和越界访问
valgrind --leak-check=full --show-leak-kinds=all ./build/openblas_demo

# 仅检查错误（更快速）
valgrind --tool=memcheck ./build/openblas_demo --ntt
```

### 使用 AddressSanitizer

在 CMakeLists.txt 的 Debug 选项中添加 ASan 支持：

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug \
         -DCMAKE_CXX_FLAGS_DEBUG="-g -O0 -fsanitize=address -fno-omit-frame-pointer"
make -j$(nproc)

# 运行（ASan 会自动检测内存问题）
./build/openblas_demo
```

### 使用 GDB 调试 OpenBLAS 内部函数

```bash
gdb ./build/openblas_demo

# 在 OpenBLAS 的 GEMM 内核上设置断点
# (gdb) break dgemm_kernel        # GEMM 微内核（函数名因 CPU 架构而异）
# (gdb) break cblas_dgemm         # CBLAS 接口层

# 查看 OpenBLAS 使用的内核名称
# (gdb) break cblas_dgemm
# (gdb) run --bench
# (gdb) step                      # 进入 cblas_dgemm
# (gdb) info sharedlibrary        # 查看加载的共享库
```

### 性能分析

```bash
# 使用 perf 进行性能分析（Linux）
perf record -g ./build/openblas_demo --bench
perf report

# 使用 gprof（需编译时添加 -pg）
cmake .. -DCMAKE_CXX_FLAGS="-pg"
make -j$(nproc)
./build/openblas_demo --bench
gprof ./build/openblas_demo gmon.out | head -50
```

---

## 代码说明

### BLAS 模块 (`include/blas_demo.h`, `src/blas_demo.cpp`)

封装了 OpenBLAS 的 CBLAS C 接口，提供类型安全的 C++ 封装：

| 函数 | BLAS 函数 | 级别 | 说明 |
|------|----------|------|------|
| `dot_product()` | `cblas_ddot` | Level 1 | 向量点积 |
| `axpy()` | `cblas_daxpy` | Level 1 | y = αx + y |
| `scal()` | `cblas_dscal` | Level 1 | x = αx |
| `vector_norm()` | `cblas_dnrm2` | Level 1 | L2 范数 |
| `vector_asum()` | `cblas_dasum` | Level 1 | L1 范数 |
| `gemv()` | `cblas_dgemv` | Level 2 | 矩阵-向量乘法 |
| `gemm()` | `cblas_dgemm` | Level 3 | 通用矩阵乘法 |
| `symm()` | `cblas_dsymm` | Level 3 | 对称矩阵乘法 |

所有矩阵使用 **列优先（Column-major）** 存储，与 Fortran/LAPACK 惯例一致。

### NTT 模块 (`include/ntt.h`, `src/ntt.cpp`)

实现了完整的数论变换，支持两组密码学标准参数：

| 参数集 | 模数 q | 维度 n | 单位根 ω | 对应方案 |
|--------|--------|--------|---------|---------|
| KyberParams | 3329 | 256 | 17 | Kyber / ML-KEM |
| DilithiumParams | 8380417 | 256 | 1753 | Dilithium / ML-DSA |

核心算法：
- **正向 NTT**: Cooley-Tukey 蝶形（迭代版本，含比特反转置换）
- **逆向 NTT**: Gentleman-Sande 蝶形 + 归一化
- **多项式乘法**: NTT → 逐点相乘 → INTT，复杂度 O(n log n)
- **朴素乘法**: O(n²) 的参考实现，用于正确性验证

### 数学背景

NTT 在格密码中用于加速多项式环 `Z_q[x]/(x^n + 1)` 上的乘法：

```
朴素乘法: O(n²) —— 对于 n=256 约 65536 次模乘法
NTT 乘法: O(n log n) —— 对于 n=256 约 2048 次模乘法 + O(n) 逐点乘法
```

---

## 常见问题

### Q: 编译时提示找不到 `cblas.h`？

A: 确认已安装 `libopenblas-dev`（Debian/Ubuntu）或 `openblas-devel`（CentOS/RHEL）。如果安装在非标准路径，设置环境变量：

```bash
export OPENBLAS_DIR=/opt/OpenBLAS
```

### Q: 运行时提示找不到 `libopenblas.so`？

A: 更新动态链接器缓存：

```bash
sudo ldconfig
# 或指定库搜索路径
export LD_LIBRARY_PATH=/opt/OpenBLAS/lib:$LD_LIBRARY_PATH
```

### Q: 如何控制 OpenBLAS 使用的线程数？

A: 通过环境变量或 C API：

```bash
# 使用 4 个线程
export OPENBLAS_NUM_THREADS=4
./build/openblas_demo

# 使用所有核心
export OPENBLAS_NUM_THREADS=$(nproc)
./build/openblas_demo
```

在代码中调用：

```cpp
openblas_set_num_threads(4);
```

### Q: 如何添加自己的测试用例？

A: 在 `src/main.cpp` 中添加新的 `demo_*` 函数，并在 `main()` 中注册即可。

### Q: NTT 实现可以直接用于生产环境吗？

A: **不建议**。本实现侧重教学可读性，生产环境请使用经过充分审计的库：
- [OpenFHE](https://github.com/openfheorg/openfhe-development)（FHE）
- [liboqs](https://github.com/open-quantum-safe/liboqs)（后量子密码）
- [PQClean](https://github.com/PQClean/PQClean)（PQC 参考实现）

---

## 许可证

本项目仅供学习和研究使用。OpenBLAS 遵循 BSD 3-Clause 许可证。
