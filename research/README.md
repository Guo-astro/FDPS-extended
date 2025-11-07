# FDPS Research Projects

This directory contains production-level research projects using the FDPS framework.

## Directory Structure

```
research/
├── 00-templates/          # Project templates
│   └── basic-sph-template/  # Standard SPH project template
├── 01-your-project/       # Your research projects (numbered)
├── 02-another-project/
└── README.md              # This file
```

## Getting Started

### 1. Create a New Project

Copy the template to start a new research project:

```bash
cd /path/to/FDPS/research
cp -r 00-templates/basic-sph-template 01-my-stellar-collision-study
cd 01-my-stellar-collision-study
```

### 2. Customize Your Project

Edit the following files:
- `README.md` - Project description and goals
- `NOTES.md` - Research log and decisions  
- `src/main.cpp` - Main simulation code
- `Makefile` - Update PROJECT name

### 3. Build and Run

```bash
# Using nix (recommended for reproducibility)
nix develop /path/to/FDPS --command bash -c 'make clean && make'

# Or directly if MPI is in PATH
make

# Run simulation
./scripts/run_simulation.sh
```

### 4. Analyze Results

```bash
# Analyze a specific run
python analysis/analyze_run.py output/run_20250107_120000

# Create animations using fdps-animator
cd ../../analysis/fdps-animator
uv run python fdps_animator.py ../../research/01-my-project/output/latest animation.mp4
```

## Project Naming Convention

Use descriptive names with optional numbering:

**Option 1: Numbered (recommended for tracking)**
```
01-stellar-collision-study/
02-galaxy-merger-simulation/
03-planet-formation-disk/
```

**Option 2: Topic-based**
```
stellar-collisions/
galaxy-mergers/
planet-formation/
```

## Best Practices

### 1. Version Control

**Keep in git:**
- Source code (`src/`)
- Scripts (`scripts/`, `analysis/`)
- Parameters (`data/parameters/`)
- Small reference data
- Documentation, notes
- Paper drafts

**Don't commit (add to .gitignore):**
- Build artifacts (`build/`, `*.o`, `*.out`)
- Large output files (`output/`)
- Temporary files

### 2. Documentation

- Maintain `README.md` with project overview
- Log decisions and observations in `NOTES.md`
- Document code with meaningful comments
- Save parameters with each run for reproducibility

### 3. Reproducibility

- Use `nix develop` for consistent environment
- Save all parameters with output
- Log git commit hash with runs
- Keep analysis scripts versioned

### 4. Testing

- Write unit tests for physics modules
- Validate against analytical solutions
- Perform convergence studies
- Document validation in `NOTES.md`

### 5. Data Management

- Organize output by date: `output/run_YYYYMMDD_HHMMSS/`
- Keep parameters with each run
- Archive important runs
- Use HDF5 for large datasets

## Standard Project Structure

Every research project should follow this structure:

```
project-name/
├── README.md              # Project overview
├── NOTES.md              # Research log
├── Makefile              # Build configuration
├── src/                  # Source code
│   ├── main.cpp
│   ├── particle_types.hpp
│   └── physics/          # Physics modules
├── scripts/              # Automation
│   └── run_simulation.sh
├── analysis/             # Post-processing
│   ├── analyze_run.py
│   └── notebooks/        # Jupyter notebooks
├── data/                 # Input data
│   ├── parameters/
│   ├── initial/
│   └── reference/
├── output/               # Simulation output (gitignored)
├── figures/              # Publication figures
├── tests/                # Project tests
└── paper/                # LaTeX drafts
```

## Using FDPS Features

### Unit Systems

The framework provides type-safe unit systems:

```cpp
using namespace PS::UnitSystem;

// Choose your unit system
Mass<GalacticUnit> m_galaxy(1.0);  // 1 solar mass
Length<GalacticUnit> r_galaxy(1.0); // 1 parsec
Velocity<GalacticUnit> v_galaxy(1.0); // 1 km/s

// Or use simulation units
Mass<SimulationUnit> m_sim(1.0);

// Convert between systems
Mass<CGSUnit> m_cgs = convertUnit<CGSUnit>(m_galaxy);
```

### Output Formats

Use the flexible output system:

```cpp
#include "fdps_output.hpp"

PS::Output::OutputMetadata metadata;
metadata.time = current_time;
metadata.total_particles = n_particles;

// CSV output
auto writer = PS::Output::createOutputWriter(
    PS::Output::OutputFormat::CSV,
    "output/snapshot",
    metadata
);

writer->open();
writer->writeParticles<ParticleType, GalacticUnit>(
    particles, n_local
);
writer->close();
```

### Visualization

Use the Python animator in `analysis/fdps-animator/`:

```python
from fdps_animator import FDPSTimeSeries, FDPSAnimator

# Load data
ts = FDPSTimeSeries()
ts.load_directory("output/run_20250107", pattern="*.csv")

# Create animation
animator = FDPSAnimator(ts)
animator.create_2d_scatter_animation(
    output_file="evolution.mp4",
    x_col='pos_x',
    y_col='pos_y',
    color_col='dens',
    fps=15
)
```

## HPC Cluster Usage

For running on HPC clusters, create a submission script:

```bash
#!/bin/bash
#SBATCH --job-name=fdps_sim
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=output/job_%j.log

module load openmpi
module load hdf5

cd $SLURM_SUBMIT_DIR
mpirun -np $SLURM_NTASKS ./research_project.out
```

## Troubleshooting

### Build fails: "mpicxx not found"

**Solution:** Enter nix environment or install MPI:
```bash
# Using nix
nix develop /path/to/FDPS

# Or install MPI
brew install open-mpi  # macOS
apt install libopenmpi-dev  # Ubuntu/Debian
```

### Linking errors with HDF5

**Solution:** The Makefile auto-detects HDF5. If not found:
```bash
# Install HDF5
brew install hdf5  # macOS
apt install libhdf5-dev  # Ubuntu/Debian

# Or use nix
nix develop /path/to/FDPS
```

## Resources

- **FDPS Documentation**: `/path/to/FDPS/docs/`
- **Unit System Guide**: `/path/to/FDPS/docs/UNIT_SYSTEM_AND_OUTPUT.md`
- **Project Structure Guide**: `/path/to/FDPS/docs/RESEARCH_PROJECT_STRUCTURE.md`
- **Sample Code**: `/path/to/FDPS/sample/c++/sph/`
- **Animator**: `/path/to/FDPS/analysis/fdps-animator/`

## Support

For questions or issues:
1. Check the documentation in `docs/`
2. Review sample code in `sample/c++/`
3. Consult `NOTES.md` in template projects
4. Check FDPS GitHub repository

## Citation

When publishing research using FDPS, please cite:
- FDPS paper: [add citation]
- Any specific methods/modules you developed

## License

All code follows the FDPS license (see root LICENSE file).
