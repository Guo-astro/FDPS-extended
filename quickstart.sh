#!/usr/bin/env bash
# Quick start script for FDPS with Nix

set -e

echo "🚀 FDPS Quick Start Script"
echo "=========================="
echo ""

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo "❌ Error: Not in a git repository"
    echo "Please run this from the FDPS repository root"
    exit 1
fi

# Check if nix is installed
if ! command -v nix &> /dev/null; then
    echo "❌ Nix is not installed"
    echo ""
    echo "To install Nix, run:"
    echo "  curl -L https://nixos.org/nix/install | sh"
    echo ""
    echo "After installation, enable flakes by adding to ~/.config/nix/nix.conf:"
    echo "  experimental-features = nix-command flakes"
    exit 1
fi

echo "✅ Nix is installed"

# Check if flake.nix is tracked by git
if ! git ls-files --error-unmatch flake.nix > /dev/null 2>&1; then
    echo ""
    echo "📝 Adding flake.nix to git..."
    git add flake.nix .envrc SETUP_GUIDE.md CHANGES.md
    git add sample/c++/sph/CMakeLists.txt sample/c++/sph/conanfile.txt
    git add sample/c++/sph/build.sh sample/c++/sph/README.md
    echo "✅ Files added to git"
fi

echo ""
echo "🔨 Building FDPS SPH sample..."
echo ""

# Enter nix develop and build
nix develop --command bash -c '
    cd sample/c++/sph
    chmod +x build.sh
    ./build.sh
    
    echo ""
    echo "✅ Build complete!"
    echo ""
    echo "To run the simulation:"
    echo "  cd sample/c++/sph/build"
    echo "  ./sph.out"
    echo ""
    echo "Or with Nix:"
    echo "  nix run"
'

echo ""
echo "🎉 Setup complete! To enter the development environment:"
echo "  nix develop"
echo ""
echo "To enable automatic environment activation:"
echo "  1. Install direnv: https://direnv.net/"
echo "  2. Run: direnv allow"
echo ""
