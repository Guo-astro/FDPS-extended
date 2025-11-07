#!/usr/bin/env bash
set -e

echo "🔨 Building FDPS SPH Sample"
echo "=========================="

# Create build directory
mkdir -p build
cd build

# Install dependencies with Conan (if needed)
echo "📦 Installing dependencies with Conan..."
conan install .. --build=missing || echo "Conan install skipped (optional)"

# Configure with CMake
echo "⚙️  Configuring with CMake..."
if [ -f conan_toolchain.cmake ]; then
    cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
else
    cmake .. -DCMAKE_BUILD_TYPE=Release
fi

# Build
echo "🔧 Building..."
cmake --build . -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo ""
echo "✅ Build complete!"
echo "Executable: $(pwd)/sph.out"
