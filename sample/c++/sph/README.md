# FDPS SPH Sample - Build Guide

This directory contains a Smoothed Particle Hydrodynamics (SPH) simulation using FDPS.

## Build Options

You can build this sample using multiple methods:

### Option 1: Nix Flake (Recommended for reproducibility)

**Prerequisites:** Install [Nix](https://nixos.org/download.html) with flakes enabled.

```bash
# Enter development environment
nix develop

# Build using the provided script
./build.sh

# Or build manually with CMake
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

**With direnv (optional):**
```bash
# Install direnv, then:
direnv allow
# Environment will auto-activate when entering directory
```

### Option 2: CMake with Conan

**Prerequisites:** CMake 3.15+, Conan 1.x or 2.x

```bash
cd sample/c++/sph

# Install dependencies (if any are added to conanfile.txt)
mkdir -p build && cd build
conan install .. --build=missing

# Configure and build
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build .

# Run
./sph.out
```

### Option 3: Traditional Makefile

**Prerequisites:** C++ compiler with C++17 support

```bash
cd sample/c++/sph

# Edit Makefile to configure compiler and flags
make

# Run
./sph.out
```

## Build Configurations

### With OpenMP (Thread Parallelization)

CMake will automatically detect OpenMP and enable it if available.

For Makefile, uncomment in the Makefile:
```makefile
CC = g++
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
```

### With MPI (Distributed Parallelization)

CMake will automatically detect MPI and enable it if available.

For Makefile, uncomment in the Makefile:
```makefile
CC = mpicxx
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
```

## Running the Simulation

```bash
# Single process
./sph.out

# With MPI (if built with MPI support)
mpirun -np 4 ./sph.out

# With Nix
nix run
```

## Adding Dependencies

### With Conan

Edit `conanfile.txt` and add dependencies:
```ini
[requires]
boost/1.82.0
eigen/3.4.0

[generators]
CMakeDeps
CMakeToolchain
```

### With Nix

Edit `../../../flake.nix` and add packages to `buildInputs`:
```nix
buildInputs = with pkgs; [
  openmpi
  libomp
  boost
  eigen
];
```

## Troubleshooting

### macOS: OpenMP not found

Install OpenMP via Homebrew:
```bash
brew install libomp
```

Or use Nix:
```bash
nix develop
```

### Variable Length Arrays Error

This occurs with strict C++ compilers. The warnings are non-critical but the error about VLA is due to the FDPS library using them. Use GCC or enable C++ extensions, or use the Nix environment which provides a compatible compiler setup.

### MPI not found

Install OpenMPI:
- **macOS:** `brew install openmpi` or `nix develop`
- **Ubuntu/Debian:** `sudo apt install libopenmpi-dev`
- **Fedora:** `sudo dnf install openmpi-devel`

## Project Structure

```
.
├── main.cpp          # SPH simulation code
├── Makefile          # Traditional build system
├── CMakeLists.txt    # CMake build configuration
├── conanfile.txt     # Conan dependencies
├── build.sh          # Convenience build script
└── README.md         # This file
```

## License

See the main FDPS LICENSE file.
