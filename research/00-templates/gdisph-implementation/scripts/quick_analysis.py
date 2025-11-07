#!/usr/bin/env python3
"""
Quick analysis of GDISPH simulation results.
Creates 1D profile plots along x-axis for shock tube comparison.
"""

import sys
from pathlib import Path
import csv

def read_csv_simple(filename):
    """Simple CSV reader without pandas"""
    data = {}
    metadata = {}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find header and data start
    header_line = None
    data_start = 0
    
    for i, line in enumerate(lines):
        if line.startswith('#'):
            # Parse metadata
            if ':' in line:
                parts = line[1:].strip().split(':', 1)
                if len(parts) == 2:
                    metadata[parts[0].strip()] = parts[1].strip()
            data_start = i + 1
        else:
            header_line = line.strip().split(',')
            data_start = i + 1
            break
    
    # Initialize data arrays
    for col in header_line:
        data[col] = []
    
    # Read data
    for line in lines[data_start:]:
        if line.strip() and not line.startswith('#'):
            values = line.strip().split(',')
            for col, val in zip(header_line, values):
                try:
                    data[col].append(float(val))
                except ValueError:
                    pass
    
    return data, metadata

def plot_1d_profiles(csv_file, output_dir="plots"):
    """
    Create 1D profile plots from SPH data by binning along x-axis.
    """
    print(f"Reading: {csv_file}")
    
    # Read data
    data, metadata = read_csv_simple(csv_file)
    
    # Get data arrays
    x = data['pos_x']
    y = data['pos_y']
    z = data['pos_z']
    
    density = data['dens']
    vx = data['vel_x']
    pressure = data['pres']
    
    # Get time from metadata
    time = float(metadata.get('time', '0.0'))
    
    # Find domain bounds
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)
    z_min, z_max = min(z), max(z)
    
    print(f"\nDomain bounds:")
    print(f"  x: [{x_min:.3f}, {x_max:.3f}]")
    print(f"  y: [{y_min:.6f}, {y_max:.6f}]")
    print(f"  z: [{z_min:.6f}, {z_max:.6f}]")
    
    # Bin ALL data along x-axis (no slab filtering for 3D shock tube)
    x_filtered = x
    dens_filtered = density
    vx_filtered = vx
    p_filtered = pressure
    
    # Create bins based on actual x range
    n_bins = 100
    bins = []
    bin_edges = [x_min + (x_max - x_min) * i / n_bins for i in range(n_bins + 1)]
    
    # Print summary
    print(f"\nTime: {time:.6f}")
    print(f"Total particles: {len(x_filtered)}")
    
    for i in range(n_bins):
        x_left = bin_edges[i]
        x_right = bin_edges[i + 1]
        x_center = (x_left + x_right) / 2
        
        # Find particles in this bin
        bin_dens = []
        bin_vel = []
        bin_pres = []
        
        for j in range(len(x_filtered)):
            if x_left <= x_filtered[j] < x_right:
                bin_dens.append(dens_filtered[j])
                bin_vel.append(vx_filtered[j])
                bin_pres.append(p_filtered[j])
        
        if len(bin_dens) > 0:
            bins.append({
                'x': x_center,
                'dens': sum(bin_dens) / len(bin_dens),
                'vel': sum(bin_vel) / len(bin_vel),
                'pres': sum(bin_pres) / len(bin_pres)
            })
    
    print(f"Bins with data: {len(bins)}")
    
    print("\nProfile (sample points):")
    print("x        density    velocity   pressure")
    print("-" * 45)
    for i in range(0, len(bins), len(bins)//10 or 1):
        b = bins[i]
        print(f"{b['x']:.3f}    {b['dens']:.4f}    {b['vel']:.4f}    {b['pres']:.4f}")
    
    # Save binned data
    Path(output_dir).mkdir(exist_ok=True)
    
    snapshot_name = Path(csv_file).stem.replace('snapshot_', '').replace('_sim', '')
    output_file = Path(output_dir) / f"profile_{snapshot_name}.csv"
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['x', 'dens', 'vel', 'pres'])
        writer.writeheader()
        writer.writerows(bins)
    
    print(f"\nSaved binned profile to: {output_file}")
    
    return bins, time

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 quick_analysis.py <csv_file>")
        print("Example: python3 quick_analysis.py result/snapshot_0060_sim.csv")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    plot_1d_profiles(csv_file)
