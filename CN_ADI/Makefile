# ==========================================
# 编译器与路径配置
# ==========================================
# 如果环境变量中没有设置 LLVM_PATH 和 LIBOMP_PATH，可以在这里指定默认值
LLVM_PATH    ?= /usr/local/llvm
LIBOMP_PATH  ?= /usr/local/opt/libomp

CC           = $(LLVM_PATH)/bin/clang
CFLAGS       = -O3 -w -Xpreprocessor -fopenmp -I$(LIBOMP_PATH)/include
LDFLAGS      = -L$(LIBOMP_PATH)/lib -lomp

# ==========================================
# 目标文件与源文件（已更新 main_adv.c）
# ==========================================
TARGET       = diff_adv_sim
SRCS         = main_adv.c diff_adv_func.c
OBJS         = $(SRCS:.c=.o)
HEADERS      = diff_adv_func.h

# ==========================================
# 编译规则
# ==========================================
.PHONY: all clean run

# 默认目标
all: $(TARGET)

# 链接目标文件生成可执行文件
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS)

# 编译 main_adv.c，依赖核心头文件
main_adv.o: main_adv.c $(HEADERS)
	$(CC) $(CFLAGS) -c main_adv.c -o main_adv.o

# 编译 diff_adv_func.c
diff_adv_func.o: diff_adv_func.c $(HEADERS)
	$(CC) $(CFLAGS) -c diff_adv_func.c -o diff_adv_func.o

# ==========================================
# 辅助命令
# ==========================================
# 运行程序
run: $(TARGET)
	time ./$(TARGET)

# 清理编译生成的文件
clean:
	rm -f $(OBJS) $(TARGET)