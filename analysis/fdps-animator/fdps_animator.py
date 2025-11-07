"""
FDPS Animator - Create animations from FDPS output files

This module provides tools to read FDPS CSV and HDF5 output files
and create animations of particle simulations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
import re

try:
    import h5py
    HDF5_AVAILABLE = True
except ImportError:
    HDF5_AVAILABLE = False


class FDPSData:
    """Container for FDPS simulation data"""
    
    def __init__(self, time: float, data: pd.DataFrame, metadata: Dict[str, Any]):
        self.time = time
        self.data = data
        self.metadata = metadata
    
    def __repr__(self):
        return f"FDPSData(time={self.time}, particles={len(self.data)}, " \
               f"unit_system={self.metadata.get('unit_system', 'Unknown')})"


class FDPSReader:
    """Read FDPS output files in CSV or HDF5 format"""
    
    @staticmethod
    def read_csv(filepath: Path) -> FDPSData:
        """
        Read FDPS CSV output file with metadata
        
        Args:
            filepath: Path to CSV file
            
        Returns:
            FDPSData object containing the data and metadata
        """
        metadata = {}
        
        # Read metadata from comment lines
        with open(filepath, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                    
                # Parse metadata
                if ':' in line:
                    key, value = line[1:].split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Convert numeric values
                    try:
                        if '.' in value or 'e' in value.lower():
                            value = float(value)
                        elif value.isdigit():
                            value = int(value)
                    except (ValueError, AttributeError):
                        pass
                    
                    metadata[key] = value
        
        # Read the CSV data
        data = pd.read_csv(filepath, comment='#')
        
        # Extract time from metadata
        time = metadata.get('Time', metadata.get('time', 0.0))
        
        return FDPSData(time=time, data=data, metadata=metadata)
    
    @staticmethod
    def read_hdf5(filepath: Path) -> FDPSData:
        """
        Read FDPS HDF5 output file
        
        Args:
            filepath: Path to HDF5 file
            
        Returns:
            FDPSData object containing the data and metadata
        """
        if not HDF5_AVAILABLE:
            raise ImportError("h5py is required to read HDF5 files. Install it with: uv add h5py")
        
        metadata = {}
        
        with h5py.File(filepath, 'r') as f:
            # Read metadata attributes
            for key in f.attrs:
                metadata[key] = f.attrs[key]
            
            # Read particle data
            particles_group = f['particles']
            data_dict = {}
            
            for dataset_name in particles_group:
                dataset = particles_group[dataset_name]
                data_dict[dataset_name] = dataset[:]
                
                # Read unit attributes if available
                if 'unit' in dataset.attrs:
                    metadata[f'{dataset_name}_unit'] = dataset.attrs['unit']
        
        # Convert to DataFrame
        data = pd.DataFrame(data_dict)
        
        # Extract time
        time = metadata.get('time', metadata.get('Time', 0.0))
        
        return FDPSData(time=time, data=data, metadata=metadata)
    
    @staticmethod
    def read_file(filepath: Path) -> FDPSData:
        """
        Automatically detect file format and read
        
        Args:
            filepath: Path to output file
            
        Returns:
            FDPSData object
        """
        filepath = Path(filepath)
        
        if filepath.suffix == '.csv':
            return FDPSReader.read_csv(filepath)
        elif filepath.suffix in ['.h5', '.hdf5']:
            return FDPSReader.read_hdf5(filepath)
        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}")


class FDPSTimeSeries:
    """
    Manage a time series of FDPS snapshots
    """
    
    def __init__(self):
        self.snapshots: List[FDPSData] = []
    
    def add_snapshot(self, snapshot: FDPSData):
        """Add a snapshot to the time series"""
        self.snapshots.append(snapshot)
        # Keep sorted by time
        self.snapshots.sort(key=lambda x: x.time)
    
    def load_directory(self, directory: Path, pattern: str = "*.csv"):
        """
        Load all files matching pattern from a directory
        
        Args:
            directory: Directory containing output files
            pattern: Glob pattern to match files (default: "*.csv")
        """
        directory = Path(directory)
        files = sorted(directory.glob(pattern))
        
        print(f"Loading {len(files)} files from {directory}")
        
        for filepath in files:
            try:
                snapshot = FDPSReader.read_file(filepath)
                self.add_snapshot(snapshot)
            except Exception as e:
                print(f"Warning: Could not read {filepath}: {e}")
        
        print(f"Loaded {len(self.snapshots)} snapshots")
        print(f"Time range: {self.get_time_range()}")
    
    def get_time_range(self) -> Tuple[float, float]:
        """Get the time range of loaded snapshots"""
        if not self.snapshots:
            return (0.0, 0.0)
        return (self.snapshots[0].time, self.snapshots[-1].time)
    
    def __len__(self):
        return len(self.snapshots)
    
    def __getitem__(self, idx):
        return self.snapshots[idx]


class FDPSAnimator:
    """
    Create animations from FDPS simulation data
    """
    
    def __init__(self, time_series: FDPSTimeSeries, figsize: Tuple[float, float] = (10, 8)):
        """
        Initialize animator
        
        Args:
            time_series: FDPSTimeSeries object with loaded snapshots
            figsize: Figure size in inches (width, height)
        """
        self.time_series = time_series
        self.figsize = figsize
        self.fig = None
        self.ax = None
    
    def create_2d_scatter_animation(
        self,
        output_file: str,
        x_col: str = 'pos_x',
        y_col: str = 'pos_y',
        color_col: Optional[str] = None,
        size_col: Optional[str] = None,
        xlim: Optional[Tuple[float, float]] = None,
        ylim: Optional[Tuple[float, float]] = None,
        fps: int = 10,
        dpi: int = 100,
        title: str = "FDPS Simulation",
        cmap: str = 'viridis',
    ):
        """
        Create a 2D scatter plot animation
        
        Args:
            output_file: Output filename (e.g., 'animation.mp4' or 'animation.gif')
            x_col: Column name for x-coordinate
            y_col: Column name for y-coordinate
            color_col: Column name for color mapping (optional)
            size_col: Column name for size mapping (optional)
            xlim: X-axis limits (auto if None)
            ylim: Y-axis limits (auto if None)
            fps: Frames per second
            dpi: Resolution
            title: Animation title
            cmap: Colormap name
        """
        if len(self.time_series) == 0:
            raise ValueError("No snapshots loaded")
        
        # Determine limits if not provided
        if xlim is None or ylim is None:
            all_x = []
            all_y = []
            for snapshot in self.time_series:
                all_x.extend(snapshot.data[x_col].values)
                all_y.extend(snapshot.data[y_col].values)
            
            if xlim is None:
                x_margin = (max(all_x) - min(all_x)) * 0.1
                xlim = (min(all_x) - x_margin, max(all_x) + x_margin)
            if ylim is None:
                y_margin = (max(all_y) - min(all_y)) * 0.1
                ylim = (min(all_y) - y_margin, max(all_y) + y_margin)
        
        # Determine color limits if color column is used
        vmin, vmax = None, None
        if color_col:
            all_colors = []
            for snapshot in self.time_series:
                all_colors.extend(snapshot.data[color_col].values)
            vmin, vmax = min(all_colors), max(all_colors)
        
        # Setup figure
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        
        # Initialize scatter plot
        snapshot0 = self.time_series[0]
        x = snapshot0.data[x_col].values
        y = snapshot0.data[y_col].values
        
        scatter_kwargs = {'alpha': 0.6}
        
        if color_col:
            colors = snapshot0.data[color_col].values
            scatter_kwargs['c'] = colors
            scatter_kwargs['vmin'] = vmin
            scatter_kwargs['vmax'] = vmax
            scatter_kwargs['cmap'] = cmap
        
        if size_col:
            sizes = snapshot0.data[size_col].values
            # Normalize sizes
            sizes = (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-10) * 100 + 10
            scatter_kwargs['s'] = sizes
        else:
            scatter_kwargs['s'] = 10
        
        scatter = self.ax.scatter(x, y, **scatter_kwargs)
        
        # Add colorbar if using colors
        if color_col:
            cbar = plt.colorbar(scatter, ax=self.ax)
            unit = snapshot0.metadata.get(f'{color_col}_unit', '')
            if unit:
                cbar.set_label(f'{color_col} [{unit}]')
            else:
                cbar.set_label(color_col)
        
        # Setup axes
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_xlabel(f'{x_col} [{snapshot0.metadata.get(f"{x_col}_unit", "")}]')
        self.ax.set_ylabel(f'{y_col} [{snapshot0.metadata.get(f"{y_col}_unit", "")}]')
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        
        # Title and time display
        time_text = self.ax.text(
            0.02, 0.98, '',
            transform=self.ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )
        self.ax.set_title(title)
        
        def update(frame_idx):
            """Update function for animation"""
            snapshot = self.time_series[frame_idx]
            
            x = snapshot.data[x_col].values
            y = snapshot.data[y_col].values
            
            # Update positions
            offsets = np.column_stack([x, y])
            scatter.set_offsets(offsets)
            
            # Update colors if applicable
            if color_col:
                colors = snapshot.data[color_col].values
                scatter.set_array(colors)
            
            # Update sizes if applicable
            if size_col:
                sizes = snapshot.data[size_col].values
                sizes = (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-10) * 100 + 10
                scatter.set_sizes(sizes)
            
            # Update time display
            time_text.set_text(f'Time: {snapshot.time:.3f}\nFrame: {frame_idx + 1}/{len(self.time_series)}')
            
            return scatter, time_text
        
        # Create animation
        anim = animation.FuncAnimation(
            self.fig,
            update,
            frames=len(self.time_series),
            interval=1000 / fps,
            blit=False
        )
        
        # Save animation
        output_path = Path(output_file)
        print(f"Saving animation to {output_path}")
        
        if output_path.suffix == '.gif':
            anim.save(output_path, writer='pillow', fps=fps, dpi=dpi)
        elif output_path.suffix in ['.mp4', '.avi', '.mov']:
            anim.save(output_path, writer='ffmpeg', fps=fps, dpi=dpi)
        else:
            raise ValueError(f"Unsupported output format: {output_path.suffix}")
        
        print(f"Animation saved successfully!")
        plt.close(self.fig)
    
    def create_histogram_animation(
        self,
        output_file: str,
        column: str,
        bins: int = 50,
        fps: int = 10,
        dpi: int = 100,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
    ):
        """
        Create an animated histogram
        
        Args:
            output_file: Output filename
            column: Column to histogram
            bins: Number of bins
            fps: Frames per second
            dpi: Resolution
            title: Plot title
            xlabel: X-axis label
        """
        if len(self.time_series) == 0:
            raise ValueError("No snapshots loaded")
        
        # Determine data range across all snapshots
        all_data = []
        for snapshot in self.time_series:
            all_data.extend(snapshot.data[column].values)
        
        data_min, data_max = min(all_data), max(all_data)
        
        # Setup figure
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        
        if title is None:
            title = f'Distribution of {column}'
        if xlabel is None:
            unit = self.time_series[0].metadata.get(f'{column}_unit', '')
            xlabel = f'{column} [{unit}]' if unit else column
        
        self.ax.set_title(title)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel('Count')
        self.ax.grid(True, alpha=0.3)
        
        time_text = self.ax.text(
            0.98, 0.98, '',
            transform=self.ax.transAxes,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )
        
        def update(frame_idx):
            """Update function for animation"""
            self.ax.clear()
            self.ax.set_title(title)
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel('Count')
            self.ax.grid(True, alpha=0.3)
            
            snapshot = self.time_series[frame_idx]
            data = snapshot.data[column].values
            
            n, bins_edges, patches = self.ax.hist(
                data, bins=bins, range=(data_min, data_max),
                alpha=0.7, edgecolor='black'
            )
            
            time_text_new = self.ax.text(
                0.98, 0.98,
                f'Time: {snapshot.time:.3f}\nFrame: {frame_idx + 1}/{len(self.time_series)}',
                transform=self.ax.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
            )
            
            return patches
        
        # Create animation
        anim = animation.FuncAnimation(
            self.fig,
            update,
            frames=len(self.time_series),
            interval=1000 / fps,
            blit=False
        )
        
        # Save animation
        output_path = Path(output_file)
        print(f"Saving animation to {output_path}")
        
        if output_path.suffix == '.gif':
            anim.save(output_path, writer='pillow', fps=fps, dpi=dpi)
        elif output_path.suffix in ['.mp4', '.avi', '.mov']:
            anim.save(output_path, writer='ffmpeg', fps=fps, dpi=dpi)
        else:
            raise ValueError(f"Unsupported output format: {output_path.suffix}")
        
        print(f"Animation saved successfully!")
        plt.close(self.fig)


def main():
    """Example usage"""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python fdps_animator.py <output_directory> [output_file]")
        print("\nExample:")
        print("  python fdps_animator.py ../result animation.mp4")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "animation.mp4"
    
    # Load time series
    time_series = FDPSTimeSeries()
    time_series.load_directory(output_dir, pattern="*.csv")
    
    if len(time_series) == 0:
        print("No snapshot files found!")
        sys.exit(1)
    
    # Create animator
    animator = FDPSAnimator(time_series)
    
    # Create 2D animation
    animator.create_2d_scatter_animation(
        output_file,
        x_col='pos_x',
        y_col='pos_y',
        color_col='dens',
        fps=10,
        title="FDPS SPH Simulation"
    )


if __name__ == '__main__':
    main()
