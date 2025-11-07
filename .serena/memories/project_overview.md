# FDPS Project Overview

## Purpose
FDPS (Framework for Developing Particle Simulators) is a general-purpose, high-performance library for particle simulations. It provides efficient algorithms and data structures for N-body simulations, including tree-based methods for force calculations, domain decomposition for parallel computing, and interfaces for C++, C, and Fortran.

Current version: 8.0a

## Tech Stack
- **Primary Language**: C++ (modern C++ with RAII, smart pointers, etc.)
- **Supported Languages**: C, Fortran (via interfaces)
- **Build System**: CMake with Ninja generator
- **Parallel Computing**: OpenMPI, OpenMP
- **Package Management**: Conan (for dependencies)
- **Development Environment**: Nix (flakes for reproducible builds)
- **Analysis Tools**: Python (with scripts in analysis/ folder)
- **Testing**: Unit tests (framework like GoogleTest implied)
- **Linting**: clang-tidy
- **Documentation**: PDF manuals in doc/ folder

## Codebase Structure
- `src/`: Core library headers and implementation
- `sample/`: Example programs in C, C++, Fortran
- `tests/`: Unit test files
- `analysis/`: Python scripts for data analysis and visualization
- `scripts/`: Utility scripts (e.g., interface generators)
- `build/`: Build artifacts (CMake, Ninja files)
- `doc/`: Documentation and papers
- `pikg/`: Related project or submodule

## Key Features
- Tree-based force calculation algorithms
- Domain decomposition for MPI parallelization
- Particle-mesh methods
- Morton key ordering
- Reallocatable arrays and memory pools
- Time profiling and timers

## Development Workflow
1. Use Nix for reproducible environment
2. Build with CMake/Ninja
3. Test with unit tests
4. Lint with clang-tidy
5. Analyze with Python scripts