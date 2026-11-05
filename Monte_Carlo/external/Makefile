CXX = c++
PYTHIA_CONFIG = ./pythia8318/bin/pythia8-config
PYTHIA_FLAGS = $(shell $(PYTHIA_CONFIG) --cxxflags --libs)

.PHONY: all run
all: pp_fixed_target

pp_fixed_target: pp_fixed_target.cpp
	$(CXX) $< $(PYTHIA_FLAGS) -o $@

pp_single_event: pp_single_event.cpp
	$(CXX) $< $(PYTHIA_FLAGS) -o $@

run: pp_fixed_target
	./pp_fixed_target
