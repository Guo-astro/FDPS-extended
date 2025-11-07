{
  description = "FDPS - Framework for Developing Particle Simulators";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        
        # Development shell with all dependencies
        devShell = pkgs.mkShell {
          buildInputs = with pkgs; [
            # Compilers and build tools
            gcc
            cmake
            gnumake
            pkg-config
            
            # Parallel computing libraries
            openmpi
            llvmPackages.openmp
            
            # Conan package manager
            conan
            
            # Additional tools
            git
            
            # Optional: for visualization and analysis
            python3
            python3Packages.numpy
            python3Packages.matplotlib
          ];

          shellHook = ''
            echo "🚀 FDPS Development Environment"
            echo "================================"
            echo "Available compilers:"
            echo "  - gcc: $(gcc --version | head -n1)"
            echo "  - g++: $(g++ --version | head -n1)"
            echo ""
            echo "Available tools:"
            echo "  - cmake: $(cmake --version | head -n1)"
            echo "  - mpicc: $(command -v mpicc && mpicc --version | head -n1 || echo 'not found')"
            echo "  - conan: $(conan --version 2>/dev/null || echo 'not found')"
            echo ""
            echo "To build the SPH sample with CMake:"
            echo "  cd sample/c++/sph"
            echo "  mkdir -p build && cd build"
            echo "  conan install .. --build=missing"
            echo "  cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release"
            echo "  cmake --build ."
            echo ""
            echo "To build with Make (traditional):"
            echo "  cd sample/c++/sph"
            echo "  make"
            echo ""
          '';

          # Set environment variables
          CC = "${pkgs.gcc}/bin/gcc";
          CXX = "${pkgs.gcc}/bin/g++";
          OMPI_CC = "${pkgs.gcc}/bin/gcc";
          OMPI_CXX = "${pkgs.gcc}/bin/g++";
        };

        # Package definition for FDPS SPH sample
        fdps-sph = pkgs.stdenv.mkDerivation {
          pname = "fdps-sph";
          version = "8.0a";
          
          src = ./.;
          
          nativeBuildInputs = with pkgs; [
            cmake
            pkg-config
          ];
          
          buildInputs = with pkgs; [
            openmpi
            llvmPackages.openmp
          ];
          
          configurePhase = ''
            cd sample/c++/sph
            mkdir -p build
            cd build
            cmake .. \
              -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_INSTALL_PREFIX=$out
          '';
          
          buildPhase = ''
            cmake --build .
          '';
          
          installPhase = ''
            mkdir -p $out/bin
            cp sph.out $out/bin/
          '';
          
          meta = with pkgs.lib; {
            description = "FDPS SPH Sample - Smoothed Particle Hydrodynamics simulation";
            homepage = "https://github.com/FDPS/FDPS";
            license = licenses.mit;
            platforms = platforms.unix;
          };
        };

      in
      {
        # Development shell
        devShells.default = devShell;
        
        # Package
        packages = {
          default = fdps-sph;
          fdps-sph = fdps-sph;
        };
        
        # Apps
        apps.default = {
          type = "app";
          program = "${fdps-sph}/bin/sph.out";
        };
      }
    );
}
