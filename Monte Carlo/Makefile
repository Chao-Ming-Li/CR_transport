# Homebrew LLVM supports OpenMP on macOS.
LLVM_PREFIX := $(shell brew --prefix llvm 2>/dev/null)
ifeq ($(origin CXX),default)
CXX := $(if $(LLVM_PREFIX),$(LLVM_PREFIX)/bin/clang++,clang++)
endif

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
OPENMP_FLAGS ?= -fopenmp
THREADS ?= $(shell sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
TARGET ?= transport
TARGETS := transport spallation hadronic_loss energy_loss transport_spallation_loss

.PHONY: all run run-transport-spallation-loss clean

all: $(TARGET)

transport: test.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(OPENMP_FLAGS) $(LDLIBS) -o $@

spallation: spallation.cpp cross_section.hpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(OPENMP_FLAGS) $(LDLIBS) -o $@

hadronic_loss: hadronic_loss.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(OPENMP_FLAGS) $(LDLIBS) -o $@

energy_loss: energy_loss.cpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(OPENMP_FLAGS) $(LDLIBS) -o $@

transport_spallation_loss: transport_spallation_loss.cpp cross_section.hpp Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $< $(LDFLAGS) $(OPENMP_FLAGS) $(LDLIBS) -o $@

run: $(TARGET)
	OMP_NUM_THREADS=$(THREADS) time ./$(TARGET)

run-transport-spallation-loss: transport_spallation_loss
	OMP_NUM_THREADS=$(THREADS) time ./transport_spallation_loss

clean:
	$(RM) $(TARGETS)
