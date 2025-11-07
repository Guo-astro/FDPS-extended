# Production Research Project Folder Structure Recommendation

## Current FDPS Structure Analysis

```
FDPS/
├── src/               # FDPS library source (framework code)
├── sample/            # Example/tutorial code
│   ├── c/
│   ├── c++/
│   │   └── sph/       # SPH sample (tutorial level)
│   └── fortran/
├── tests/             # Unit tests for FDPS framework
├── analysis/          # Analysis tools
│   ├── fdps-animator/ # Python visualization package
│   ├── plots/
│   ├── scripts/
│   └── animations/
├── docs/              # Framework documentation
├── build/             # Build artifacts (gitignored)
└── papers/            # Research papers
```

## Recommended Production Research Structure

### Option 1: Separate Top-Level Research Directory (Recommended)

```
FDPS/
├── src/               # FDPS library (DO NOT MODIFY)
├── sample/            # Examples (DO NOT MODIFY)
├── tests/             # Framework tests
├── analysis/          # Shared analysis tools
├── docs/              # Framework docs
│
├── research/          # 🎯 NEW: Your production research projects
│   │
│   ├── 00-templates/  # Project templates and boilerplate
│   │   ├── basic-sph/
│   │   ├── nbody/
│   │   └── README.md
│   │
│   ├── 01-stellar-collision/  # Example project 1
│   │   ├── src/
│   │   │   ├── main.cpp
│   │   │   ├── initial_conditions.hpp
│   │   │   ├── physics/
│   │   │   │   ├── eos.hpp
│   │   │   │   └── gravity.hpp
│   │   │   └── io/
│   │   │       ├── output_manager.hpp
│   │   │       └── checkpoint.hpp
│   │   ├── scripts/
│   │   │   ├── run_suite.sh
│   │   │   ├── postprocess.py
│   │   │   └── convergence_test.py
│   │   ├── analysis/
│   │   │   ├── visualize.py
│   │   │   ├── energy_evolution.py
│   │   │   └── notebooks/
│   │   │       ├── 01_initial_conditions.ipynb
│   │   │       ├── 02_evolution_analysis.ipynb
│   │   │       └── 03_final_state.ipynb
│   │   ├── data/
│   │   │   ├── initial/       # Initial condition files
│   │   │   ├── parameters/    # Parameter files
│   │   │   └── reference/     # Reference data for validation
│   │   ├── output/            # Simulation output (gitignored)
│   │   ├── figures/           # Publication-ready figures
│   │   ├── build/             # Build artifacts (gitignored)
│   │   ├── tests/             # Project-specific tests
│   │   │   ├── unit/
│   │   │   ├── integration/
│   │   │   └── convergence/
│   │   ├── CMakeLists.txt
│   │   ├── Makefile
│   │   ├── README.md          # Project documentation
│   │   ├── NOTES.md           # Research notes
│   │   └── paper/             # LaTeX paper draft
│   │       ├── main.tex
│   │       ├── figures/
│   │       └── bibliography.bib
│   │
│   ├── 02-galaxy-merger/      # Example project 2
│   │   └── [same structure as above]
│   │
│   ├── 03-planet-formation/   # Example project 3
│   │   └── [same structure as above]
│   │
│   └── README.md              # Research projects index
│
└── .gitignore                 # Updated to ignore research/*/output/, research/*/build/
```

### Option 2: External Research Repository (For Multiple Researchers)

```
fdps-research-projects/        # Separate git repository
├── README.md
├── .gitmodules                # FDPS as submodule
├── FDPS/                      # Git submodule → your FDPS repo
│
├── john-stellar/              # Researcher 1's projects
│   ├── collision-study/
│   └── binary-evolution/
│
├── jane-galactic/             # Researcher 2's projects
│   ├── merger-simulations/
│   └── dark-matter-halos/
│
└── shared/                    # Shared utilities
    ├── analysis-tools/
    ├── plotting-styles/
    └── validation-data/
```

## Best Practices for Production Research

### 1. Project Naming Convention

```
research/
├── 01-descriptive-name/    # Number prefix for ordering
├── 02-another-project/
└── 03-latest-work/
```

Or by topic:
```
research/
├── stellar-collisions/
├── galaxy-mergers/
└── planet-formation/
```

### 2. Standard Project Structure

Every research project should have:

```
project-name/
├── README.md              # Purpose, how to build, how to run
├── NOTES.md              # Research log, decisions, observations
├── CMakeLists.txt        # Build configuration
├── Makefile              # Alternative build (optional)
├── src/                  # Source code
│   ├── main.cpp
│   ├── particle_types.hpp
│   ├── initial_conditions.hpp
│   ├── physics/          # Physics modules
│   ├── io/               # I/O modules
│   └── utils/            # Utilities
├── scripts/              # Automation scripts
│   ├── run_simulation.sh
│   ├── submit_job.sh     # HPC cluster submission
│   └── parameter_sweep.py
├── analysis/             # Post-processing
│   ├── plot_results.py
│   ├── calculate_stats.py
│   └── notebooks/        # Jupyter notebooks
├── data/                 # Input data
│   ├── initial/
│   ├── parameters/
│   └── reference/
├── output/               # Simulation output (gitignored)
├── figures/              # Publication figures
├── tests/                # Tests
│   ├── unit_tests/
│   └── validation/
└── paper/                # Paper drafts
    ├── main.tex
    └── figures/
```

### 3. Recommended CMakeLists.txt Template

```cmake
cmake_minimum_required(VERSION 3.15)
project(stellar_collision LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FDPS include path
set(FDPS_ROOT "${CMAKE_SOURCE_DIR}/../../.." CACHE PATH "FDPS root directory")
include_directories(${FDPS_ROOT}/src)

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS})

# Find OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Optional: HDF5
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    add_definitions(-DPARTICLE_SIMULATOR_USE_HDF5)
endif()

# Compilation flags
add_compile_options(
    -O3
    -ffast-math
    -funroll-loops
    -DPARTICLE_SIMULATOR_THREAD_PARALLEL
    -DPARTICLE_SIMULATOR_MPI_PARALLEL
)

# Source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Executable
add_executable(simulation ${SOURCES})

target_link_libraries(simulation
    MPI::MPI_CXX
    $<$<BOOL:${HDF5_FOUND}>:${HDF5_CXX_LIBRARIES}>
    $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_CXX>
)

# Install
install(TARGETS simulation DESTINATION bin)
```

### 4. Recommended Makefile Template

```makefile
# Project Configuration
PROJECT = stellar_collision
FDPS_ROOT = ../../..
PS_PATH = -I $(FDPS_ROOT)/src/

# Compiler Detection
MPICXX := $(shell command -v mpicxx 2> /dev/null)
ifndef MPICXX
    $(error mpicxx not found! Please install MPI or use: nix develop)
endif

CC = mpicxx
CFLAGS = -std=c++17 -O3 -ffast-math -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# Optional HDF5
HDF5_FOUND := $(shell pkg-config --exists hdf5 && echo yes)
ifeq ($(HDF5_FOUND),yes)
    CFLAGS += -DPARTICLE_SIMULATOR_USE_HDF5
    LIBS += $(shell pkg-config --libs hdf5)
endif

# Source files
SOURCES = $(wildcard src/*.cpp src/*/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
PROGRAM = $(PROJECT).out

.PHONY: all clean run test

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	@echo "Linking..."
	@$(CC) $(CFLAGS) $(OBJECTS) -o $(PROGRAM) $(LIBS)
	@echo "✓ Build complete: $(PROGRAM)"

%.o: %.cpp
	@echo "Compiling $<..."
	@$(CC) -c $< $(CFLAGS) $(PS_PATH) -o $@

clean:
	@rm -f $(OBJECTS) $(PROGRAM)
	@echo "✓ Cleaned"

run: $(PROGRAM)
	@mkdir -p output
	@mpirun -np 4 ./$(PROGRAM)

test: $(PROGRAM)
	@echo "Running tests..."
	@$(MAKE) -C tests run
```

### 5. Git Strategy

#### research/.gitignore additions:

```gitignore
# Build artifacts
*/build/
*/output/
*/*.out
*/*.o

# Data files (large)
*/data/initial/*.dat
*/data/initial/*.h5
*/output/**/*.h5
*/output/**/*.csv

# Temporary files
*/.ipynb_checkpoints/
*/__pycache__/
```

#### Keep in git:
- Source code (`src/`)
- Scripts (`scripts/`, `analysis/`)
- Parameters (`data/parameters/`)
- Small reference data
- Documentation
- Paper drafts

#### Don't commit:
- Build artifacts
- Large simulation output
- Temporary files

### 6. Using with Nix

Add to your `flake.nix`:

```nix
packages = {
  # ... existing packages ...

  research-stellar = pkgs.stdenv.mkDerivation {
    name = "fdps-stellar-collision";
    src = ./research/01-stellar-collision;
    
    buildInputs = with pkgs; [
      gcc
      openmpi
      hdf5
      llvmPackages.openmp
    ];
    
    buildPhase = ''
      make -j$NIX_BUILD_CORES
    '';
    
    installPhase = ''
      mkdir -p $out/bin
      cp *.out $out/bin/
    '';
  };
};
```

### 7. Running Simulations

Create a standard run script:

```bash
#!/bin/bash
# scripts/run_simulation.sh

set -e

# Configuration
NPROCS=8
OUTPUT_DIR="output/run_$(date +%Y%m%d_%H%M%S)"
PARAM_FILE="data/parameters/default.param"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Copy parameters for reproducibility
cp "$PARAM_FILE" "$OUTPUT_DIR/parameters.param"

# Run simulation
mpirun -np $NPROCS ./stellar_collision.out \
    --param "$PARAM_FILE" \
    --output "$OUTPUT_DIR" \
    | tee "$OUTPUT_DIR/simulation.log"

echo "✓ Simulation complete: $OUTPUT_DIR"
```

### 8. Analysis Pipeline

```python
# analysis/analyze_run.py
#!/usr/bin/env python3
"""
Analyze simulation output and generate figures
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent / 'analysis' / 'fdps-animator'))

from fdps_animator import FDPSTimeSeries, FDPSAnimator
import pandas as pd
import matplotlib.pyplot as plt

def analyze_run(output_dir):
    """Analyze a single simulation run"""
    output_path = Path(output_dir)
    
    # Load time series
    ts = FDPSTimeSeries()
    ts.load_directory(output_path, pattern="*.csv")
    
    # Create animations
    animator = FDPSAnimator(ts)
    animator.create_2d_scatter_animation(
        output_file=str(output_path / "evolution.mp4"),
        x_col='pos_x',
        y_col='pos_y',
        color_col='dens',
        fps=15
    )
    
    # Energy evolution
    energies = []
    for snap in ts:
        total_energy = snap.data['eng'].sum()
        energies.append((snap.time, total_energy))
    
    df = pd.DataFrame(energies, columns=['time', 'energy'])
    
    plt.figure(figsize=(10, 6))
    plt.plot(df['time'], df['energy'])
    plt.xlabel('Time')
    plt.ylabel('Total Energy')
    plt.title('Energy Conservation')
    plt.grid(True)
    plt.savefig(output_path / 'energy_evolution.png', dpi=300)
    
    print(f"✓ Analysis complete: {output_dir}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python analyze_run.py <output_directory>")
        sys.exit(1)
    
    analyze_run(sys.argv[1])
```

## Summary & Recommendation

**For your production research, I recommend:**

1. **Create `/Users/guo-opt-p148/FDPS/research/` directory** for all production work
2. **Keep `sample/` untouched** - it's for learning and examples
3. **Structure each project** with the template above
4. **Use consistent naming** (numbered or topic-based)
5. **Leverage the shared `analysis/` tools** (fdps-animator, etc.)
6. **Version control** your research projects separately if they get large
7. **Use nix develop** for reproducible builds

This structure:
- ✅ Separates production from samples
- ✅ Keeps each project self-contained
- ✅ Makes it easy to share/publish individual projects
- ✅ Reuses FDPS library without modification
- ✅ Follows scientific computing best practices
- ✅ Works well with HPC clusters
- ✅ Makes reproducibility straightforward

Would you like me to create the initial `research/` directory structure with a template project?
