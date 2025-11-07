PS_PATH = -I ../../../src/

# Compiler selection
CC = mpicxx
CFLAGS = -std=c++17 -O3 -ffast-math -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

SRC_DIR = src
OUTPUT_DIR = result
ANIM_DIR = animations
ANALYSIS_DIR = analysis
CPPOBJS = $(patsubst $(SRC_DIR)/%.cpp, %.o, $(wildcard $(SRC_DIR)/*.cpp))
CPPHDRS = $(wildcard $(SRC_DIR)/*.h)
PROGRAM = sph.out
NIX_FLAKE = ../../..
ANIMATOR_PATH = ../../../analysis/fdps-animator

# MPI check - only enforced for targets that need it
MPICXX := $(shell command -v mpicxx 2> /dev/null)
NEEDS_MPI = all %.o run

.PHONY:	clean all run simulate animate compare-analytic distclean help check-mpi

# MPI availability check target
check-mpi:
ifndef MPICXX
	$(error mpicxx not found! Please install OpenMPI or MPICH, or use 'make build-nix' to build with Nix environment)
endif

all:	check-mpi $(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	$(SRC_DIR)/%.cpp $(CPPHDRS) check-mpi
	@echo "Building $< ..."
	@$(CC) -c $< $(CFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

# Build with nix environment (recommended)
build-nix:
	@echo "🔨 Building with Nix environment..."
	@nix develop $(NIX_FLAKE) --command bash -c "make all"

# Run simulation with nix environment
run-nix:
	@echo "🚀 Running simulation with Nix environment..."
	@mkdir -p $(OUTPUT_DIR)
	@nix develop $(NIX_FLAKE) --command bash -c "./$(PROGRAM)"

# Complete workflow: build, run, and generate animations
simulate: build-nix run-nix animate
	@echo "✨ Complete simulation workflow finished!"
	@echo "📊 Results in: $(OUTPUT_DIR)/"
	@echo "🎬 Animations in: $(ANIM_DIR)/"

# Generate animations from simulation output
animate:
	@echo "🎬 Generating animations..."
	@mkdir -p $(ANIM_DIR)
	@cd $(ANIMATOR_PATH) && uv run python fdps_animator.py \
		$(CURDIR)/$(OUTPUT_DIR) \
		$(CURDIR)/$(ANIM_DIR)/sph_density.gif
	@echo "✅ Animation created: $(ANIM_DIR)/sph_density.gif"

# Run analytic comparison
compare-analytic:
	@echo "📈 Comparing with analytic solution..."
	@mkdir -p $(ANALYSIS_DIR)
	@cd $(ANIMATOR_PATH) && uv run python $(CURDIR)/scripts/compare_analytic.py
	@echo "✅ Comparison plots created in $(ANALYSIS_DIR)/"

run:	all
	@echo "Running simulation..."
	@mkdir -p $(OUTPUT_DIR)
	@./$(PROGRAM)

clean:
	-rm -f *.out *.o

distclean: clean
	-rm -rf $(OUTPUT_DIR) $(ANIM_DIR) $(ANALYSIS_DIR)

# Help target
help:
	@echo "FDPS SPH Simulation Makefile"
	@echo "============================"
	@echo ""
	@echo "⚠️  This simulation requires MPI (mpicxx compiler)"
	@echo "   Recommended: Use 'make simulate' or 'make build-nix' with Nix environment"
	@echo "   Alternative: Install MPI locally (brew install open-mpi)"
	@echo ""
	@echo "Main targets:"
	@echo "  make simulate          - Complete workflow: build, run, and animate (with Nix)"
	@echo "  make build-nix         - Build using Nix development environment"
	@echo "  make run-nix           - Run simulation using Nix environment"
	@echo "  make animate           - Generate animations from output"
	@echo "  make compare-analytic  - Compare simulation with analytic solution"
	@echo ""
	@echo "Basic targets (require local MPI installation):"
	@echo "  make all              - Build the program (requires mpicxx)"
	@echo "  make run              - Run the simulation"
	@echo "  make clean            - Remove compiled files"
	@echo "  make distclean        - Remove all generated files"
	@echo ""
	@echo "Output directories:"
	@echo "  $(OUTPUT_DIR)/       - Simulation results (CSV files)"
	@echo "  $(ANIM_DIR)/         - Animation files (GIF/MP4)"
	@echo "  $(ANALYSIS_DIR)/     - Analysis and comparison plots"

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
