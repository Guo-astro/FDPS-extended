# Suggested Commands for FDPS Development

## Environment Setup
- `nix develop`: Enter the Nix development environment with all dependencies
- `git add flake.nix .envrc`: Track Nix files (required for Nix flakes)

## Building
- `cd sample/c++/sph && ./build.sh`: Build the SPH sample (convenience script)
- `mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build .`: Standard CMake build
- `ninja`: Build using Ninja (if configured)

## Testing
- `ctest`: Run tests (if CMake test targets are configured)
- Unit tests are in `tests/` folder, build and run manually or via CMake

## Linting and Code Quality
- `clang-tidy <file>`: Lint individual files (configure with project's clang-tidy profile)
- Treat warnings as errors: compile with `-Werror -Wall -Wextra -Wpedantic`

## Analysis and Visualization
- `cd analysis && python read_fdps_csv.py`: Run Python analysis scripts
- Various plotting and animation scripts in `analysis/`

## Utility Commands (macOS/Darwin)
- `ls -la`: List files
- `find . -name "*.hpp"`: Find header files
- `grep -r "pattern" src/`: Search in source code
- `git status`: Check git status
- `git diff`: See changes
- `git add . && git commit -m "message"`: Commit changes

## Running Samples
- `cd sample/c++/sph && ./build/sph.out`: Run SPH simulation
- Similar for other samples in `sample/c/`, `sample/fortran/`

## Documentation
- PDFs in `doc/` folder: tutorials, specifications
- `README.md`, `SETUP_GUIDE.md` for setup
- `CHANGES.md` for changelog

## Cleaning
- Remove build artifacts: `rm -rf build/`
- Clean temporary files: `.bak`, `.orig`, etc. before committing