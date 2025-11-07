# FDPS with Conan and Nix - Setup Guide

This guide explains how to build and run FDPS (Framework for Developing Particle Simulators) using modern build tools.

## 📦 What's Been Added

### 1. **Nix Flake** (`flake.nix`)
A declarative, reproducible development environment with:
- GCC compiler
- CMake build system
- OpenMPI for parallel computing
- Conan package manager
- Python tools for analysis

### 2. **Conan Configuration** (`sample/c++/sph/conanfile.txt`)
Package manager configuration for C++ dependencies.

### 3. **CMake Build System** (`sample/c++/sph/CMakeLists.txt`)
Modern build configuration with:
- Automatic OpenMP detection
- Automatic MPI detection
- Cross-platform support

### 4. **Build Script** (`sample/c++/sph/build.sh`)
Convenience script for quick builds.

## 🚀 Quick Start

### Option 1: Using Nix (Recommended)

**Requirements:** Nix with flakes enabled

```bash
# 1. Add flake.nix to git (Nix requires tracked files)
git add flake.nix .envrc

# 2. Enter development environment
nix develop

# 3. Build the SPH sample
cd sample/c++/sph
./build.sh

# 4. Run the simulation
./build/sph.out
```

### Option 2: Using CMake directly

**Requirements:** CMake 3.15+, C++17 compiler

```bash
cd sample/c++/sph
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./sph.out
```

### Option 3: Using Conan + CMake

**Requirements:** Conan, CMake, C++17 compiler

```bash
cd sample/c++/sph
mkdir -p build && cd build

# Install dependencies (if added to conanfile.txt)
conan install .. --build=missing

# Build
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
cmake --build .
./sph.out
```

## 🔧 Installation Prerequisites

### macOS

```bash
# Install Nix (recommended)
curl -L https://nixos.org/nix/install | sh

# Enable flakes (add to ~/.config/nix/nix.conf or /etc/nix/nix.conf)
experimental-features = nix-command flakes

# Or install dependencies via Homebrew
brew install cmake gcc openmpi conan
```

### Linux (Ubuntu/Debian)

```bash
# Install Nix (recommended)
curl -L https://nixos.org/nix/install | sh

# Or install dependencies via apt
sudo apt update
sudo apt install build-essential cmake libopenmpi-dev python3-pip
pip3 install conan
```

### Linux (Fedora)

```bash
# Install Nix (recommended)
curl -L https://nixos.org/nix/install | sh

# Or install dependencies via dnf
sudo dnf install gcc-c++ cmake openmpi-devel python3-pip
pip3 install conan
```

## 🛠️ Known Issues and Workarounds

### Issue 1: FDPS Library Compilation Error

The FDPS library has a non-const reference binding issue with Apple Clang:

```
error: non-const lvalue reference to type 'Vector3<...>' cannot bind to a temporary
```

**Workaround:** Use GCC instead of Clang:

```bash
# With Nix (automatically uses GCC)
nix develop

# With Homebrew on macOS
brew install gcc
export CC=/usr/local/bin/gcc-13
export CXX=/usr/local/bin/g++-13
```

### Issue 2: OpenMP Not Found on macOS

**Solution with Nix:**
```bash
nix develop  # OpenMP is included
```

**Solution with Homebrew:**
```bash
brew install libomp
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"
```

### Issue 3: MPI Not Found

**Solution with Nix:**
```bash
nix develop  # OpenMPI is included
```

**Solution with package managers:**
- **macOS:** `brew install openmpi`
- **Ubuntu:** `sudo apt install libopenmpi-dev`
- **Fedora:** `sudo dnf install openmpi-devel`

## 📚 Using Nix Features

### Development Shell

```bash
# Enter dev environment
nix develop

# With direnv (auto-activate on cd)
echo "use flake" > .envrc
direnv allow
```

### Building the Package

```bash
# Build the FDPS SPH package
nix build .#fdps-sph

# Run directly
nix run
```

### Updating Dependencies

```bash
# Update all inputs to latest versions
nix flake update

# Update specific input
nix flake lock --update-input nixpkgs
```

## 📦 Adding Dependencies

### With Conan

Edit `sample/c++/sph/conanfile.txt`:

```ini
[requires]
boost/1.82.0
fmt/10.1.0
spdlog/1.12.0

[generators]
CMakeDeps
CMakeToolchain

[options]
boost:shared=False
```

Then rebuild:
```bash
cd sample/c++/sph/build
conan install .. --build=missing
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake
cmake --build .
```

### With Nix

Edit `flake.nix` and add to `buildInputs`:

```nix
buildInputs = with pkgs; [
  openmpi
  libomp
  boost
  fmt
  spdlog
];
```

Then:
```bash
nix develop
```

### With CMake

Edit `sample/c++/sph/CMakeLists.txt`:

```cmake
find_package(Boost REQUIRED)
target_link_libraries(sph.out Boost::boost)
```

## 🎯 Parallel Execution

### With OpenMP (Thread Parallelization)

The build system automatically detects and enables OpenMP if available.

```bash
# Set number of threads
export OMP_NUM_THREADS=4
./sph.out
```

### With MPI (Distributed Parallelization)

```bash
# Run on 4 processes
mpirun -np 4 ./sph.out

# With Nix
nix develop -c mpirun -np 4 ./build/sph.out
```

## 🧪 Testing the Setup

```bash
# 1. Enter Nix environment
nix develop

# 2. Verify compilers
gcc --version
g++ --version
cmake --version

# 3. Build SPH sample
cd sample/c++/sph
./build.sh

# 4. Run simulation
./build/sph.out
```

## 📁 Directory Structure

```
FDPS/
├── flake.nix                 # Nix flake configuration
├── .envrc                    # direnv configuration
├── src/                      # FDPS library headers
└── sample/c++/sph/
    ├── main.cpp              # SPH simulation
    ├── Makefile              # Traditional build
    ├── CMakeLists.txt        # CMake configuration
    ├── conanfile.txt         # Conan dependencies
    ├── build.sh              # Build script
    └── README.md             # Sample documentation
```

## 🔗 Resources

- **FDPS Documentation:** Check `doc/` directory
- **Nix Flakes:** https://nixos.wiki/wiki/Flakes
- **Conan:** https://conan.io/
- **CMake:** https://cmake.org/

## 💡 Tips

1. **Use Nix for reproducibility:** Guaranteed to work across different systems
2. **Use direnv:** Automatically activates environment when entering directory
3. **Use build.sh:** Simplifies the build process
4. **Version control:** Remember to `git add flake.nix` before using Nix commands

## 🤝 Contributing

When adding new dependencies:
1. Add to `conanfile.txt` for Conan users
2. Add to `flake.nix` for Nix users
3. Update `CMakeLists.txt` to find and link the package
4. Update this documentation

## 📝 License

FDPS is provided under the MIT License. See LICENSE file for details.
