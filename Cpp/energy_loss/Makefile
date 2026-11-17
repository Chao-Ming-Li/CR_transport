CXX = c++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra -Wpedantic

ifeq ($(shell uname -s),Darwin)
LIBOMP_PREFIX ?= $(shell if command -v brew >/dev/null 2>&1; then brew --prefix libomp; fi)
OMP_CXXFLAGS ?= -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
OMP_LDLIBS ?= -L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib -lomp
else
OMP_CXXFLAGS ?= -fopenmp
OMP_LDLIBS ?= -fopenmp
endif

SOURCES = configuration.cpp field.cpp initialization.cpp spatial_solver.cpp momentum_solver.cpp nuclear_solver.cpp
HEADERS = configuration.hpp field.hpp initialization.hpp spatial_solver.hpp momentum_solver.hpp nuclear_solver.hpp cross_section.hpp
TESTS = test_field test_spatial test_momentum test_diffusion_model test_nuclear

.PHONY: all openmp test test-openmp clean
all: cr_solver
openmp: cr_solver_omp

cr_solver: main.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) main.cpp $(SOURCES) -o $@

$(TESTS): %: %.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $< $(SOURCES) -o $@

cr_solver_omp: main.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OMP_CXXFLAGS) main.cpp $(SOURCES) $(LDFLAGS) $(OMP_LDLIBS) -o $@

$(addsuffix _omp,$(TESTS)): %_omp: %.cpp $(SOURCES) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OMP_CXXFLAGS) $< $(SOURCES) $(LDFLAGS) $(OMP_LDLIBS) -o $@

test-openmp: $(addsuffix _omp,$(TESTS))
	set -e; for test in $^; do OMP_NUM_THREADS=$${OMP_NUM_THREADS:-4} ./$$test; done

test: $(TESTS)
	./test_field
	./test_spatial
	./test_momentum
	./test_diffusion_model
	./test_nuclear

clean:
	rm -f cr_solver cr_solver_omp $(TESTS) $(addsuffix _omp,$(TESTS))
