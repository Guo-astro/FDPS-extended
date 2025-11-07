#!/usr/bin/env python3
"""
Example Python script to read and analyze FDPS CSV output files
with unit metadata support.

This demonstrates how the CSV format can be easily read and processed
by analysis tools while preserving unit information.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re


class FDPSCSVReader:
    """
    Reader for FDPS CSV output files with unit metadata.
    """
    
    def __init__(self, filename):
        self.filename = filename
        self.metadata = {}
        self.units = {}
        self.data = None
        self._read_file()
    
    def _read_file(self):
        """Read the CSV file and extract metadata, units, and data."""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        # Parse metadata from comment lines
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('#'):
                if ':' in line:
                    # Extract metadata
                    key_value = line[1:].strip().split(':', 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        self.metadata[key.strip()] = value.strip()
                
                # Check for units line
                if 'Units:' in line:
                    unit_str = line.split('Units:')[1].strip()
                    self.units = dict(zip(
                        lines[i-1].strip().split(','),
                        unit_str.split(',')
                    ))
                
                data_start = i + 1
            else:
                # First non-comment line is column headers
                if not self.data:
                    break
        
        # Read the actual data
        self.data = pd.read_csv(self.filename, comment='#')
    
    def get_metadata(self, key):
        """Get a metadata value by key."""
        return self.metadata.get(key)
    
    def get_unit(self, column):
        """Get the unit for a specific column."""
        return self.units.get(column, 'unknown')
    
    def get_column(self, column):
        """Get data for a specific column as numpy array."""
        return self.data[column].values
    
    def get_position(self):
        """Get position vectors as (N, 3) array."""
        return np.column_stack([
            self.get_column('pos_x'),
            self.get_column('pos_y'),
            self.get_column('pos_z')
        ])
    
    def get_velocity(self):
        """Get velocity vectors as (N, 3) array."""
        return np.column_stack([
            self.get_column('vel_x'),
            self.get_column('vel_y'),
            self.get_column('vel_z')
        ])
    
    def summary(self):
        """Print a summary of the data."""
        print(f"File: {self.filename}")
        print(f"{'='*60}")
        print("\nMetadata:")
        for key, value in self.metadata.items():
            print(f"  {key}: {value}")
        
        print(f"\nColumns and Units:")
        for col in self.data.columns:
            unit = self.get_unit(col)
            print(f"  {col:15s} [{unit}]")
        
        print(f"\nData Shape: {self.data.shape}")
        print(f"Number of particles: {len(self.data)}")
        
        print(f"\nData Statistics:")
        print(self.data.describe())
    
    def plot_2d(self, x_col='pos_x', y_col='pos_y', color_col='dens', 
                output_file=None):
        """
        Create a 2D scatter plot of the particles.
        
        Args:
            x_col: Column name for x-axis (default: 'pos_x')
            y_col: Column name for y-axis (default: 'pos_y')
            color_col: Column name for color mapping (default: 'dens')
            output_file: If provided, save plot to this file
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        x = self.get_column(x_col)
        y = self.get_column(y_col)
        c = self.get_column(color_col)
        
        scatter = ax.scatter(x, y, c=c, s=10, alpha=0.6, cmap='viridis')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(f"{color_col} [{self.get_unit(color_col)}]")
        
        # Labels
        ax.set_xlabel(f"{x_col} [{self.get_unit(x_col)}]")
        ax.set_ylabel(f"{y_col} [{self.get_unit(y_col)}]")
        
        # Title with metadata
        time = self.get_metadata('Time')
        unit_system = self.get_metadata('Unit System')
        title = f"Particle Distribution (t={time}, {unit_system} units)"
        ax.set_title(title)
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def plot_histogram(self, column='dens', bins=50, output_file=None):
        """
        Create a histogram of a specified column.
        
        Args:
            column: Column name to plot
            bins: Number of bins
            output_file: If provided, save plot to this file
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        data = self.get_column(column)
        
        ax.hist(data, bins=bins, alpha=0.7, edgecolor='black')
        
        ax.set_xlabel(f"{column} [{self.get_unit(column)}]")
        ax.set_ylabel("Count")
        
        time = self.get_metadata('Time')
        unit_system = self.get_metadata('Unit System')
        title = f"{column} Distribution (t={time}, {unit_system} units)"
        ax.set_title(title)
        
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add statistics text
        stats_text = (
            f"Mean: {np.mean(data):.3e}\n"
            f"Std: {np.std(data):.3e}\n"
            f"Min: {np.min(data):.3e}\n"
            f"Max: {np.max(data):.3e}"
        )
        ax.text(0.98, 0.98, stats_text,
                transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {output_file}")
        else:
            plt.show()
        
        plt.close()


def compare_unit_systems(file_sim, file_cgs, file_gal):
    """
    Compare the same snapshot in different unit systems.
    
    Args:
        file_sim: Path to simulation unit file
        file_cgs: Path to CGS unit file
        file_gal: Path to galactic unit file
    """
    print("Comparing different unit systems for the same snapshot:")
    print("="*60)
    
    # Read all files
    readers = {
        'Simulation': FDPSCSVReader(file_sim),
        'CGS': FDPSCSVReader(file_cgs),
        'Galactic': FDPSCSVReader(file_gal)
    }
    
    # Check that particle IDs match
    ids_sim = readers['Simulation'].get_column('id')
    ids_cgs = readers['CGS'].get_column('id')
    ids_gal = readers['Galactic'].get_column('id')
    
    if not (np.array_equal(ids_sim, ids_cgs) and np.array_equal(ids_sim, ids_gal)):
        print("WARNING: Particle IDs don't match across files!")
        return
    
    print(f"\n✓ All files contain the same {len(ids_sim)} particles")
    
    # Compare a specific quantity (mass) in different units
    print("\nMass of first particle in different unit systems:")
    for name, reader in readers.items():
        mass = reader.get_column('mass')[0]
        unit = reader.get_unit('mass')
        print(f"  {name:12s}: {mass:.6e} [{unit}]")
    
    # Compare positions
    print("\nPosition of first particle in different unit systems:")
    for name, reader in readers.items():
        pos = reader.get_position()[0]
        unit = reader.get_unit('pos_x')
        print(f"  {name:12s}: ({pos[0]:.6e}, {pos[1]:.6e}, {pos[2]:.6e}) [{unit}]")


def main():
    """Main demonstration function."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage:")
        print(f"  {sys.argv[0]} <csv_file>                 # Read and summarize")
        print(f"  {sys.argv[0]} <csv_file> --plot          # Create visualizations")
        print(f"  {sys.argv[0]} <sim> <cgs> <gal> --compare # Compare unit systems")
        return
    
    if '--compare' in sys.argv:
        if len(sys.argv) < 5:
            print("Error: Need three files for comparison")
            return
        compare_unit_systems(sys.argv[1], sys.argv[2], sys.argv[3])
        return
    
    # Read the file
    reader = FDPSCSVReader(sys.argv[1])
    
    # Print summary
    reader.summary()
    
    # Create plots if requested
    if '--plot' in sys.argv or '-p' in sys.argv:
        print("\n" + "="*60)
        print("Creating visualizations...")
        
        # 2D position plot
        base = Path(sys.argv[1]).stem
        reader.plot_2d(output_file=f"{base}_pos2d.png")
        
        # Density histogram
        reader.plot_histogram('dens', output_file=f"{base}_dens_hist.png")
        
        # Velocity histogram
        reader.plot_histogram('vel_x', output_file=f"{base}_velx_hist.png")
        
        print("\nVisualization complete!")


if __name__ == '__main__':
    main()
