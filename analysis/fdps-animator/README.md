# FDPS Animator

Create animations from FDPS (Framework for Developing Particle Simulators) output files.

## Features

- **Multiple File Formats**: Read CSV and HDF5 output files
- **Automatic Metadata Parsing**: Extract simulation time, unit systems, and other metadata
- **Flexible Animations**: Create 2D scatter plots, histograms, and custom visualizations
- **Color & Size Mapping**: Map particle properties to colors and sizes
- **Video Export**: Export to MP4, AVI, MOV, or GIF formats
- **Modern Package Management**: Uses `uv` for fast, reliable dependency management

## Installation

This project uses [uv](https://github.com/astral-sh/uv) as the package manager.

```bash
# Navigate to the animator directory
cd /path/to/FDPS/analysis/fdps-animator

# Install dependencies (uv automatically manages the virtual environment)
uv sync
```

## Quick Start

### Basic Usage

```python
from fdps_animator import FDPSTimeSeries, FDPSAnimator

# Load all CSV files from a directory
time_series = FDPSTimeSeries()
time_series.load_directory("../../sample/c++/sph/result", pattern="*.csv")

# Create animator
animator = FDPSAnimator(time_series)

# Generate 2D scatter animation
animator.create_2d_scatter_animation(
    output_file="simulation.mp4",
    x_col='pos_x',
    y_col='pos_y',
    color_col='dens',  # Color by density
    fps=10,
    title="SPH Simulation"
)
```

### Command Line Usage

```bash
# Run from the analysis/fdps-animator directory
uv run python fdps_animator.py ../../sample/c++/sph/result animation.mp4
```

## API Reference

### FDPSData

Container for a single simulation snapshot.

**Attributes:**
- `time`: Simulation time
- `data`: pandas DataFrame with particle data
- `metadata`: Dictionary with simulation metadata

### FDPSReader

Read FDPS output files.

**Methods:**
- `read_csv(filepath)`: Read CSV file
- `read_hdf5(filepath)`: Read HDF5 file  
- `read_file(filepath)`: Auto-detect format and read

### FDPSTimeSeries

Manage a time series of snapshots.

**Methods:**
- `add_snapshot(snapshot)`: Add a single snapshot
- `load_directory(directory, pattern)`: Load all matching files from directory
- `get_time_range()`: Get (min_time, max_time) tuple
- `__len__()`: Number of snapshots
- `__getitem__(idx)`: Access snapshot by index

### FDPSAnimator

Create animations from time series data.

#### Constructor

```python
animator = FDPSAnimator(time_series, figsize=(10, 8))
```

#### Methods

##### `create_2d_scatter_animation()`

Create a 2D particle scatter plot animation.

```python
animator.create_2d_scatter_animation(
    output_file="animation.mp4",
    x_col='pos_x',              # X-coordinate column
    y_col='pos_y',              # Y-coordinate column
    color_col='dens',           # Optional: color mapping
    size_col='mass',            # Optional: size mapping
    xlim=(-10, 10),             # Optional: x-axis limits
    ylim=(-10, 10),             # Optional: y-axis limits
    fps=10,                     # Frames per second
    dpi=100,                    # Resolution
    title="Simulation",         # Plot title
    cmap='viridis'             # Colormap
)
```

##### `create_histogram_animation()`

Create an animated histogram.

```python
animator.create_histogram_animation(
    output_file="histogram.mp4",
    column='dens',              # Column to histogram
    bins=50,                    # Number of bins
    fps=10,                     # Frames per second
    dpi=100,                    # Resolution
    title="Density Distribution",
    xlabel="Density [g/cm^3]"
)
```

## Advanced Examples

### Custom Color Mapping

```python
# Color by energy, size by mass
animator.create_2d_scatter_animation(
    output_file="energy_colored.mp4",
    x_col='pos_x',
    y_col='pos_y',
    color_col='eng',
    size_col='mass',
    cmap='plasma',
    fps=15
)
```

### 3D Projection to 2D

```python
# Create XY, XZ, and YZ projections
for x, y, name in [('pos_x', 'pos_y', 'xy'), 
                    ('pos_x', 'pos_z', 'xz'),
                    ('pos_y', 'pos_z', 'yz')]:
    animator.create_2d_scatter_animation(
        output_file=f"projection_{name}.mp4",
        x_col=x,
        y_col=y,
        color_col='dens',
        title=f"{name.upper()} Projection"
    )
```

### Multiple Histograms

```python
# Create histograms for different properties
for prop in ['dens', 'pres', 'eng']:
    animator.create_histogram_animation(
        output_file=f"{prop}_histogram.mp4",
        column=prop,
        bins=50,
        fps=10
    )
```

## File Format Requirements

### CSV Format

CSV files should include:
- Metadata as comment lines starting with `#`
- Column headers
- Particle data

Example:
```
# Time: 0.5
# Total Particles: 1000
# Unit System: CGS
# Description: SPH simulation
id,mass,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,dens,eng,pres
0,1.0,0.5,0.3,0.2,0.01,0.02,0.03,1.0,2.5,1.5
...
```

### HDF5 Format

HDF5 files should include:
- Root-level attributes for metadata (time, unit_system, etc.)
- A `particles` group containing datasets for each property
- Optional `unit` attributes on datasets

## Video Export Formats

Supported formats:
- **MP4** (`.mp4`): H.264 video (requires ffmpeg)
- **AVI** (`.avi`): Uncompressed or compressed video (requires ffmpeg)
- **MOV** (`.mov`): QuickTime video (requires ffmpeg)
- **GIF** (`.gif`): Animated GIF (uses pillow, no ffmpeg required)

**Note**: For MP4/AVI/MOV export, you need `ffmpeg` installed:

```bash
# macOS
brew install ffmpeg

# Linux
sudo apt-get install ffmpeg  # Debian/Ubuntu
sudo yum install ffmpeg      # RHEL/CentOS

# Or use conda/mamba
conda install -c conda-forge ffmpeg
```

## Dependencies

Managed automatically by `uv`:

- `numpy`: Numerical computations
- `pandas`: Data manipulation
- `matplotlib`: Plotting and animation
- `h5py`: HDF5 file support
- `pillow`: GIF export

## Development

### Adding Dependencies

```bash
uv add package-name
```

### Project Structure

```
fdps-animator/
├── fdps_animator.py    # Main module
├── main.py             # Entry point
├── pyproject.toml      # Project configuration
├── README.md           # This file
└── .venv/              # Virtual environment (auto-managed by uv)
```

## Troubleshooting

### "No module named 'h5py'"

If you're reading HDF5 files and get this error, install h5py:
```bash
uv add h5py
```

### "ffmpeg not found"

For MP4/AVI/MOV export, install ffmpeg (see Video Export Formats section above).

### "No snapshots loaded"

Check that:
1. The directory path is correct
2. Files match the pattern (default: `*.csv`)
3. Files are in the correct format

### Memory Issues with Large Simulations

For very large particle counts, consider:
- Processing fewer snapshots
- Reducing the DPI
- Using GIF format instead of video

## License

This project is part of FDPS and follows the same license terms.
