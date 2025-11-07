#!/usr/bin/env python3
"""
Generate animation comparing GDISPH simulation results with Sod shock tube analytical solution.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
import glob

# Sod shock tube analytical solution
def sod_shock_analytical(x, t, gamma=1.4):
    """
    Compute analytical solution for Sod shock tube problem at time t.
    
    Initial conditions:
    - Left state (x <= 0): rho=1.0, P=1.0, v=0.0
    - Right state (x > 0): rho=0.125, P=0.1, v=0.0
    """
    # Left and right states
    rho_L, P_L, v_L = 1.0, 1.0, 0.0
    rho_R, P_R, v_R = 0.125, 0.1, 0.0
    
    # Sound speeds
    c_L = np.sqrt(gamma * P_L / rho_L)
    c_R = np.sqrt(gamma * P_R / rho_R)
    
    # Solve for star region pressure (iterative method simplified)
    # For Sod problem, analytical values are well-known
    P_star = 0.30313  # Star region pressure
    v_star = 0.92745  # Star region velocity
    
    # Wave speeds
    # Head of rarefaction
    x_head = -c_L * t
    
    # Tail of rarefaction  
    c_star_L = c_L * (P_star / P_L) ** ((gamma - 1) / (2 * gamma))
    x_tail = (v_star - c_star_L) * t
    
    # Contact discontinuity
    x_contact = v_star * t
    
    # Shock front
    v_shock = v_R + c_R * np.sqrt((gamma + 1) / (2 * gamma) * P_star / P_R + (gamma - 1) / (2 * gamma))
    x_shock = v_shock * t
    
    # Compute solution at each x location
    rho = np.zeros_like(x)
    P = np.zeros_like(x)
    v = np.zeros_like(x)
    u = np.zeros_like(x)  # Specific internal energy
    
    for i, xi in enumerate(x):
        if xi < x_head:
            # Left state (undisturbed)
            rho[i] = rho_L
            P[i] = P_L
            v[i] = v_L
        elif xi < x_tail:
            # Rarefaction fan
            c = (2 / (gamma + 1)) * (c_L + xi / t)
            rho[i] = rho_L * (c / c_L) ** (2 / (gamma - 1))
            P[i] = P_L * (c / c_L) ** (2 * gamma / (gamma - 1))
            v[i] = (2 / (gamma + 1)) * (c_L + xi / t)
        elif xi < x_contact:
            # Left star region
            rho[i] = rho_L * (P_star / P_L) ** (1 / gamma)
            P[i] = P_star
            v[i] = v_star
        elif xi < x_shock:
            # Right star region
            rho[i] = rho_R * ((P_star / P_R) + (gamma - 1) / (gamma + 1)) / ((gamma - 1) / (gamma + 1) * (P_star / P_R) + 1)
            P[i] = P_star
            v[i] = v_star
        else:
            # Right state (undisturbed)
            rho[i] = rho_R
            P[i] = P_R
            v[i] = v_R
    
    # Compute specific internal energy: u = P / ((gamma - 1) * rho)
    u = P / ((gamma - 1) * rho)
    
    return rho, P, v, u


def load_snapshot(filename):
    """Load a simulation snapshot CSV file."""
    # Skip comment lines that start with #
    df = pd.read_csv(filename, comment='#')
    return df


def create_animation(result_dir, output_file='gdisph_vs_analytical.mp4'):
    """
    Create animation comparing GDISPH simulation with analytical solution.
    """
    result_path = Path(result_dir)
    
    # Find all snapshot files
    snapshot_files = sorted(glob.glob(str(result_path / 'snapshot_*_sim.csv')))
    
    if not snapshot_files:
        print(f"No snapshot files found in {result_dir}")
        return
    
    print(f"Found {len(snapshot_files)} snapshot files")
    
    # Setup figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('GDISPH vs Analytical Solution (Sod Shock Tube)', fontsize=16)
    
    axes = axes.flatten()
    
    # Analytical solution x-grid
    x_analytical = np.linspace(-1.0, 1.0, 1000)
    
    def animate(frame_num):
        # Clear all axes
        for ax in axes:
            ax.clear()
        
        if frame_num >= len(snapshot_files):
            return
        
        # Load simulation data
        df = load_snapshot(snapshot_files[frame_num])
        
        # Extract time from filename or data
        # Assuming time is stored or can be computed
        # For now, estimate from frame number
        t = frame_num * 0.02  # Approximate time step
        
        # Compute analytical solution
        rho_analytical, P_analytical, v_analytical, u_analytical = sod_shock_analytical(x_analytical, t if t > 0 else 0.001)
        
        # Plot density
        axes[0].plot(x_analytical, rho_analytical, 'r-', label='Analytical', linewidth=2)
        axes[0].plot(df['pos_x'], df['dens'], 'b.', label='GDISPH', markersize=2, alpha=0.6)
        axes[0].set_xlabel('Position (x)')
        axes[0].set_ylabel('Density')
        axes[0].set_title('Density')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        axes[0].set_xlim(-1, 1)
        axes[0].set_ylim(0, 1.2)
        
        # Plot pressure
        axes[1].plot(x_analytical, P_analytical, 'r-', label='Analytical', linewidth=2)
        axes[1].plot(df['pos_x'], df['pres'], 'b.', label='GDISPH', markersize=2, alpha=0.6)
        axes[1].set_xlabel('Position (x)')
        axes[1].set_ylabel('Pressure')
        axes[1].set_title('Pressure')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        axes[1].set_xlim(-1, 1)
        axes[1].set_ylim(0, 1.2)
        
        # Plot velocity
        axes[2].plot(x_analytical, v_analytical, 'r-', label='Analytical', linewidth=2)
        axes[2].plot(df['pos_x'], df['vel_x'], 'b.', label='GDISPH', markersize=2, alpha=0.6)
        axes[2].set_xlabel('Position (x)')
        axes[2].set_ylabel('Velocity')
        axes[2].set_title('Velocity')
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)
        axes[2].set_xlim(-1, 1)
        axes[2].set_ylim(-0.2, 1.2)
        
        # Plot internal energy
        axes[3].plot(x_analytical, u_analytical, 'r-', label='Analytical', linewidth=2)
        axes[3].plot(df['pos_x'], df['eng'], 'b.', label='GDISPH', markersize=2, alpha=0.6)
        axes[3].set_xlabel('Position (x)')
        axes[3].set_ylabel('Specific Internal Energy')
        axes[3].set_title('Internal Energy')
        axes[3].legend()
        axes[3].grid(True, alpha=0.3)
        axes[3].set_xlim(-1, 1)
        axes[3].set_ylim(0, 3.0)
        
        # Add time annotation
        fig.text(0.5, 0.95, f'Time: {t:.4f}', ha='center', fontsize=14)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Create animation
    print("Creating animation...")
    anim = animation.FuncAnimation(fig, animate, frames=len(snapshot_files), 
                                   interval=200, repeat=True)
    
    # Save animation
    print(f"Saving animation to {output_file}...")
    anim.save(output_file, writer='ffmpeg', fps=5, dpi=100)
    print(f"Animation saved successfully!")
    
    plt.close()


if __name__ == '__main__':
    import sys
    
    result_dir = sys.argv[1] if len(sys.argv) > 1 else 'result'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'gdisph_vs_analytical.mp4'
    
    print(f"Processing simulation data from: {result_dir}")
    print(f"Output file: {output_file}")
    
    create_animation(result_dir, output_file)
