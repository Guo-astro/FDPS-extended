# FDPS Research Project Template

This is a template for creating production-level research projects with FDPS.

## Quick Start

1. Copy this template to create a new project:
   ```bash
   cd /path/to/FDPS/research
   cp -r 00-templates/basic-sph-template 01-my-research-project
   cd 01-my-research-project
   ```

2. Update `README.md` with your project details

3. Build and run:
   ```bash
   # Using nix (recommended)
   nix develop /path/to/FDPS --command bash -c 'make clean && make'
   
   # Or directly if MPI is installed
   make
   ```

4. Run simulation:
   ```bash
   ./scripts/run_simulation.sh
   ```

5. Analyze results:
   ```bash
   python analysis/analyze_run.py output/latest
   ```

## Project Structure

- `src/` - Source code
  - `main.cpp` - Main simulation code
  - `particle_types.hpp` - Particle structures
  - `initial_conditions.hpp` - IC generation
  - `physics/` - Physics modules (EOS, gravity, etc.)
  - `io/` - I/O modules (output, checkpointing)

- `scripts/` - Automation scripts
  - `run_simulation.sh` - Standard run script
  - `parameter_sweep.py` - Parameter space exploration

- `analysis/` - Post-processing
  - `analyze_run.py` - Standard analysis pipeline
  - `notebooks/` - Jupyter notebooks for exploration

- `data/` - Input data
  - `parameters/` - Parameter files
  - `initial/` - Initial condition files
  - `reference/` - Reference/validation data

- `output/` - Simulation output (gitignored)
- `figures/` - Publication-ready figures
- `tests/` - Project-specific tests
- `paper/` - LaTeX paper draft

## Build Options

### Using Makefile (simple)
```bash
make          # Build
make clean    # Clean
make run      # Build and run with 4 MPI ranks
make test     # Run tests
```

### Using CMake (flexible)
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

### Using Nix (reproducible)
```bash
nix develop /path/to/FDPS --command bash -c 'make'
```

## Configuration

Edit `data/parameters/default.param` to configure simulation parameters.

## Best Practices

1. **Version control**: Keep source code, scripts, and small data in git
2. **Documentation**: Update NOTES.md with research decisions and observations
3. **Reproducibility**: Always save parameters with output
4. **Testing**: Write tests for physics modules
5. **Checkpointing**: Save checkpoints for long runs

## Citation

If you use this research in publications, please cite:
- FDPS: [FDPS paper citation]
- Your specific modifications/methods

## License

Follow the same license as FDPS.
