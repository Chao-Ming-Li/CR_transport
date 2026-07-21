# ============================================
# Makefile for CR Diffusion Simulation
# 使用 Homebrew LLVM (clang++) + libomp
# ============================================

# Homebrew LLVM 和 libomp 路径（对应你环境变量里的设置）
LLVM_PATH := /opt/homebrew/opt/llvm
LIBOMP_PATH := /opt/homebrew/opt/libomp

# 编译器：直接指定 Homebrew LLVM 的 clang++（而不是依赖 PATH 里的顺序）
CXX := $(LLVM_PATH)/bin/clang++

# 编译标志：C++20 + 优化 + OpenMP（通过 -Xpreprocessor 让 clang 正确处理 #pragma omp）
CXXFLAGS := -std=c++20 -O3 -Wall -Wextra \
            -Xpreprocessor -fopenmp \
            -I$(LIBOMP_PATH)/include

# 链接标志：链接 libomp 动态库
LDFLAGS := -L$(LIBOMP_PATH)/lib -lomp

# 目标可执行文件名
TARGET := diffusion_sim

# 自动查找所有 .cpp 源文件
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(wildcard *.hpp)

# ============================================
# 规则
# ============================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)
	@echo "Build complete: ./$(TARGET)"

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
	@echo "Cleaned build artifacts."

run: $(TARGET)
	time ./$(TARGET) &

rebuild: clean all

.PHONY: all clean run rebuild