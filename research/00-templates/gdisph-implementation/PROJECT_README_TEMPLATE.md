# Research Project: [Your Title]

## Overview

Brief description of the research question and approach.

## Author

- Name
- Institution
- Email
- Date started: YYYY-MM-DD

## Scientific Goals

1. Investigate ...
2. Measure ...
3. Compare ...

## Methods

- SPH simulation with FDPS
- Physical model: [e.g., self-gravity + hydrodynamics]
- Equation of state: [e.g., ideal gas]
- Resolution: [e.g., 10^5 particles]

## Parameters

See `data/parameters/default.param` for simulation parameters.

Key parameters:
- N particles: 100,000
- Box size: 10 kpc
- Timestep: adaptive
- Duration: 1 Gyr

## Build Instructions

### Prerequisites

- MPI (OpenMPI or MPICH)
- OpenMP
- Optional: HDF5

### Quick Start

```bash
# Using Nix (recommended for reproducibility)
cd /path/to/FDPS
nix develop

cd research/your-project-name
make
./scripts/run_simulation.sh
```

### Without Nix

```bash
# Ensure mpicxx is in PATH
make
```

## Running Simulations

### Single Run

```bash
./scripts/run_simulation.sh
```

### Parameter Sweep

```bash
python scripts/parameter_sweep.py
```

### HPC Cluster

```bash
sbatch scripts/submit_job.sh
```

## Analysis

```bash
# Analyze a run
python analysis/analyze_run.py output/run_20250107_120000

# Create animations
python analysis/create_animations.py output/run_20250107_120000

# Interactive exploration
jupyter lab analysis/notebooks/
```

## Results

Key findings:
1. ...
2. ...

## Publications

- Paper 1: [Citation]
- Conference: [Citation]

## Data

- Raw simulation data: [location]
- Processed data: [location]
- Figures: `figures/`

## Changelog

### 2025-01-07
- Initial project setup
- Implemented basic SPH

### 2025-01-XX
- Added ...

## References

1. Reference 1
2. Reference 2

## License

Same as FDPS (see LICENSE in root directory)
