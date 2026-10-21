# Use Homebrew LLVM for OpenMP on macOS; use g++ on Linux.
ifeq ($(shell uname -s),Darwin)
LLVM_PREFIX := $(shell brew --prefix llvm 2>/dev/null)
ifeq ($(origin CXX),default)
CXX := $(if $(LLVM_PREFIX),$(LLVM_PREFIX)/bin/clang++,clang++)
endif
else
ifeq ($(origin CXX),default)
CXX := g++
endif
endif

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
OPENMP_FLAGS ?= -fopenmp
THREADS ?= $(shell nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
TARGET := diffusion

.PHONY: all run clean

all: $(TARGET)

$(TARGET): diffusion.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(OPENMP_FLAGS) $(LDLIBS) -o $@

run: $(TARGET)
	OMP_NUM_THREADS=$(THREADS) time ./$(TARGET) &

clean:
	$(RM) $(TARGET)
