# ============================================================================
# OpenBLAS 演示项目 - CMake 构建系统
# ============================================================================
# 本项目演示 BLAS 基本操作（Level 1/2/3）以及数论变换（NTT）
# 适用于密码学研究方向（格密码、同态加密等）
# ============================================================================

cmake_minimum_required(VERSION 3.14)
project(OpenBLAS_Demo
    VERSION 1.0.0
    DESCRIPTION "OpenBLAS 演示项目：BLAS 基础算子与 NTT 实现"
    LANGUAGES CXX
)

# ----------------------------------------------------------------------------
# 编译标准与选项
# ----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 默认构建类型为 Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "构建类型" FORCE)
endif()

# Debug 模式下启用调试符号和 AddressSanitizer 选项
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# ----------------------------------------------------------------------------
# 查找 OpenBLAS 库
# ----------------------------------------------------------------------------
# 方式一：通过 pkg-config 查找（推荐）
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(OPENBLAS QUIET openblas)
endif()

# 方式二：手动搜索 OpenBLAS
if(NOT OPENBLAS_FOUND)
    # 常见的 OpenBLAS 安装路径
    set(OPENBLAS_SEARCH_PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib/x86_64-linux-gnu
        /opt/OpenBLAS/lib
        $ENV{OPENBLAS_DIR}/lib
    )

    find_library(OPENBLAS_LIBRARIES
        NAMES openblas blas
        PATHS ${OPENBLAS_SEARCH_PATHS}
    )

    # 查找头文件
    set(OPENBLAS_INCLUDE_SEARCH_PATHS
        /usr/include
        /usr/local/include
        /opt/OpenBLAS/include
        $ENV{OPENBLAS_DIR}/include
    )

    find_path(OPENBLAS_INCLUDE_DIRS
        NAMES cblas.h openblas/cblas.h
        PATHS ${OPENBLAS_INCLUDE_SEARCH_PATHS}
    )

    if(OPENBLAS_LIBRARIES AND OPENBLAS_INCLUDE_DIRS)
        set(OPENBLAS_FOUND TRUE)
    endif()
endif()

if(NOT OPENBLAS_FOUND)
    message(FATAL_ERROR
        "未找到 OpenBLAS 库！请安装 OpenBLAS 或设置 OPENBLAS_DIR 环境变量。\n"
        "安装方法：\n"
        "  Ubuntu/Debian: sudo apt-get install libopenblas-dev\n"
        "  CentOS/RHEL:   sudo yum install openblas-devel\n"
        "  macOS:         brew install openblas\n"
        "  源码编译:      git clone https://github.com/OpenMathLib/OpenBLAS && make && sudo make install"
    )
endif()

message(STATUS "OpenBLAS 库: ${OPENBLAS_LIBRARIES}")
message(STATUS "OpenBLAS 头文件: ${OPENBLAS_INCLUDE_DIRS}")

# ----------------------------------------------------------------------------
# 头文件包含路径
# ----------------------------------------------------------------------------
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${OPENBLAS_INCLUDE_DIRS}
)

# ----------------------------------------------------------------------------
# 库文件搜索路径
# ----------------------------------------------------------------------------
link_directories(
    ${CMAKE_SOURCE_DIR}/lib
)

# ----------------------------------------------------------------------------
# 构建目标
# ----------------------------------------------------------------------------

# 1. BLAS 演示库（静态库）
add_library(blas_demo STATIC
    src/blas_demo.cpp
)
target_include_directories(blas_demo PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${OPENBLAS_INCLUDE_DIRS}
)

# 2. NTT 库（静态库）
add_library(ntt STATIC
    src/ntt.cpp
)
target_include_directories(ntt PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# 3. 主程序可执行文件
add_executable(openblas_demo
    src/main.cpp
)
target_link_libraries(openblas_demo
    blas_demo
    ntt
    ${OPENBLAS_LIBRARIES}
    pthread
    m
)

# ----------------------------------------------------------------------------
# 安装规则
# ----------------------------------------------------------------------------
install(TARGETS openblas_demo DESTINATION bin)
install(TARGETS blas_demo ntt DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

# ----------------------------------------------------------------------------
# 辅助目标：清理构建目录
# ----------------------------------------------------------------------------
add_custom_target(clean_build
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMENT "清理构建目录"
)
