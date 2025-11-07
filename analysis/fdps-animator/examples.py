"""
Example: Create animations from FDPS output files

This script demonstrates various ways to create animations
from FDPS simulation output data.
"""

from fdps_animator import FDPSTimeSeries, FDPSAnimator
from pathlib import Path


def example_basic_animation():
    """Basic 2D scatter plot animation"""
    print("=" * 60)
    print("Example 1: Basic 2D Scatter Animation")
    print("=" * 60)
    
    # Load time series
    time_series = FDPSTimeSeries()
    time_series.load_directory("../../sample/c++/sph/result", pattern="*.csv")
    
    if len(time_series) == 0:
        print("No data files found. Run an SPH simulation first!")
        return
    
    # Create animator
    animator = FDPSAnimator(time_series, figsize=(10, 8))
    
    # Create basic animation
    animator.create_2d_scatter_animation(
        output_file="basic_animation.mp4",
        x_col='pos_x',
        y_col='pos_y',
        fps=10,
        title="Basic SPH Simulation"
    )
    
    print("✓ Created: basic_animation.mp4\n")


def example_colored_animation():
    """Animation with color mapping"""
    print("=" * 60)
    print("Example 2: Color-Mapped Animation")
    print("=" * 60)
    
    time_series = FDPSTimeSeries()
    time_series.load_directory("../../sample/c++/sph/result", pattern="*.csv")
    
    if len(time_series) == 0:
        print("No data files found!")
        return
    
    animator = FDPSAnimator(time_series, figsize=(12, 10))
    
    # Color by density
    animator.create_2d_scatter_animation(
        output_file="density_colored.mp4",
        x_col='pos_x',
        y_col='pos_y',
        color_col='dens',
        cmap='plasma',
        fps=10,
        title="SPH Simulation - Density Colored"
    )
    
    print("✓ Created: density_colored.mp4\n")


def example_histogram():
    """Animated histogram"""
    print("=" * 60)
    print("Example 3: Animated Histogram")
    print("=" * 60)
    
    time_series = FDPSTimeSeries()
    time_series.load_directory("../../sample/c++/sph/result", pattern="*.csv")
    
    if len(time_series) == 0:
        print("No data files found!")
        return
    
    animator = FDPSAnimator(time_series, figsize=(10, 6))
    
    # Density histogram
    animator.create_histogram_animation(
        output_file="density_histogram.mp4",
        column='dens',
        bins=50,
        fps=10,
        title="Density Distribution Over Time"
    )
    
    print("✓ Created: density_histogram.mp4\n")


def example_multiple_projections():
    """Create multiple 2D projections of 3D data"""
    print("=" * 60)
    print("Example 4: Multiple 2D Projections")
    print("=" * 60)
    
    time_series = FDPSTimeSeries()
    time_series.load_directory("../../sample/c++/sph/result", pattern="*.csv")
    
    if len(time_series) == 0:
        print("No data files found!")
        return
    
    animator = FDPSAnimator(time_series, figsize=(10, 10))
    
    # Create XY, XZ, and YZ projections
    projections = [
        ('pos_x', 'pos_y', 'XY'),
        ('pos_x', 'pos_z', 'XZ'),
        ('pos_y', 'pos_z', 'YZ')
    ]
    
    for x_col, y_col, name in projections:
        animator.create_2d_scatter_animation(
            output_file=f"projection_{name.lower()}.mp4",
            x_col=x_col,
            y_col=y_col,
            color_col='dens',
            fps=10,
            title=f"{name} Projection - Density Colored"
        )
        print(f"✓ Created: projection_{name.lower()}.mp4")
    
    print()


def example_gif_export():
    """Export as GIF instead of video"""
    print("=" * 60)
    print("Example 5: GIF Export (no ffmpeg required)")
    print("=" * 60)
    
    time_series = FDPSTimeSeries()
    time_series.load_directory("../../sample/c++/sph/result", pattern="*.csv")
    
    if len(time_series) == 0:
        print("No data files found!")
        return
    
    animator = FDPSAnimator(time_series, figsize=(8, 6))
    
    # Export as GIF
    animator.create_2d_scatter_animation(
        output_file="simulation.gif",
        x_col='pos_x',
        y_col='pos_y',
        color_col='dens',
        fps=5,  # Lower FPS for smaller file size
        dpi=80,  # Lower DPI for smaller file size
        title="SPH Simulation"
    )
    
    print("✓ Created: simulation.gif\n")


def main():
    """Run all examples"""
    print("\n")
    print("╔" + "=" * 58 + "╗")
    print("║" + " " * 12 + "FDPS Animator - Examples" + " " * 22 + "║")
    print("╚" + "=" * 58 + "╝")
    print()
    
    # Check if output directory exists
    result_dir = Path("../../sample/c++/sph/result")
    if not result_dir.exists():
        print("⚠ Output directory not found!")
        print(f"   Expected: {result_dir.absolute()}")
        print("   Please run an SPH simulation first to generate output files.")
        print()
        return
    
    # Run examples
    examples = [
        ("1", "Basic Animation", example_basic_animation),
        ("2", "Color-Mapped Animation", example_colored_animation),
        ("3", "Animated Histogram", example_histogram),
        ("4", "Multiple Projections", example_multiple_projections),
        ("5", "GIF Export", example_gif_export),
    ]
    
    print("Available examples:")
    for num, name, _ in examples:
        print(f"  [{num}] {name}")
    print("  [a] Run all examples")
    print()
    
    choice = input("Select example (1-5, a, or q to quit): ").strip().lower()
    
    if choice == 'q':
        return
    elif choice == 'a':
        for num, name, func in examples:
            try:
                func()
            except Exception as e:
                print(f"✗ Error in {name}: {e}\n")
    elif choice in ['1', '2', '3', '4', '5']:
        idx = int(choice) - 1
        try:
            examples[idx][2]()
        except Exception as e:
            print(f"✗ Error: {e}\n")
    else:
        print("Invalid choice!")
    
    print("Done! Check the current directory for output files.")


if __name__ == '__main__':
    main()
