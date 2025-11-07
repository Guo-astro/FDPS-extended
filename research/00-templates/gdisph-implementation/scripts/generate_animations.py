#!/usr/bin/env python3
"""
Generate animations from FDPS SPH simulation output

Run from the basic-sph-template directory:
    cd /Users/guo-opt-p148/FDPS/analysis/fdps-animator
    uv run python ../../research/00-templates/basic-sph-template/scripts/generate_animations.py
"""

from pathlib import Path
from fdps_animator import FDPSTimeSeries, FDPSAnimator

def main():
    # Paths relative to fdps-animator directory
    result_dir = Path("../../research/00-templates/basic-sph-template/result")
    anim_dir = Path("../../research/00-templates/basic-sph-template/animations")
    anim_dir.mkdir(exist_ok=True)
    
    print("=" * 60)
    print("FDPS SPH Animation Generator")
    print("=" * 60)
    
    # Load time series
    print(f"\n📁 Loading snapshots from: {result_dir}")
    time_series = FDPSTimeSeries()
    time_series.load_directory(result_dir, pattern="snapshot_*_cgs.csv")
    
    if len(time_series) == 0:
        print("❌ No snapshots found!")
        return 1
    
    print(f"✅ Loaded {len(time_series)} snapshots")
    print(f"⏱️  Time range: {time_series.get_time_range()}")
    
    # Create animator
    animator = FDPSAnimator(time_series, figsize=(12, 10))
    
    # 1. Basic XY projection with density coloring
    print("\n🎬 Creating XY projection animation...")
    animator.create_2d_scatter_animation(
        output_file=str(anim_dir / "sph_xy_density.gif"),
        x_col='pos_x',
        y_col='pos_y',
        color_col='dens',
        cmap='plasma',
        fps=5,
        dpi=100,
        title="SPH Simulation - XY Projection (Density)"
    )
    print("✅ Created: sph_xy_density.gif")
    
    # 2. XZ projection
    print("\n🎬 Creating XZ projection animation...")
    animator.create_2d_scatter_animation(
        output_file=str(anim_dir / "sph_xz_density.gif"),
        x_col='pos_x',
        y_col='pos_z',
        color_col='dens',
        cmap='plasma',
        fps=5,
        dpi=100,
        title="SPH Simulation - XZ Projection (Density)"
    )
    print("✅ Created: sph_xz_density.gif")
    
    # 3. Density histogram
    print("\n🎬 Creating density histogram animation...")
    animator.create_histogram_animation(
        output_file=str(anim_dir / "density_histogram.gif"),
        column='dens',
        bins=50,
        fps=5,
        dpi=100,
        title="Density Distribution Evolution"
    )
    print("✅ Created: density_histogram.gif")
    
    # 4. Energy colored animation
    print("\n🎬 Creating energy-colored animation...")
    animator.create_2d_scatter_animation(
        output_file=str(anim_dir / "sph_xy_energy.gif"),
        x_col='pos_x',
        y_col='pos_y',
        color_col='eng',
        cmap='viridis',
        fps=5,
        dpi=100,
        title="SPH Simulation - Energy Distribution"
    )
    print("✅ Created: sph_xy_energy.gif")
    
    print("\n" + "=" * 60)
    print("✨ All animations created successfully!")
    print(f"📂 Output directory: {anim_dir.absolute()}")
    print("=" * 60)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main())
