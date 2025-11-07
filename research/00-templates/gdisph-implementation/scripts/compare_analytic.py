#!/usr/bin/env python3
"""
Compare FDPS SPH simulation with analytic solution (Sod Shock Tube)

This script compares the numerical SPH simulation results with the
analytic Riemann solution for the 1D Sod shock tube problem.

Run from fdps-animator directory:
    cd /path/to/FDPS/analysis/fdps-animator
    uv run python /path/to/research/.../scripts/compare_analytic.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Tuple, Dict, List
import sys

# Sod shock tube initial conditions (dimensionless units)
class SodShockTube:
    """Analytic solution for 1D Sod shock tube problem"""
    
    def __init__(self):
        # Left state (x < 0.5)
        self.rho_L = 1.0
        self.p_L = 1.0
        self.u_L = 0.0
        
        # Right state (x > 0.5)
        self.rho_R = 0.125
        self.p_R = 0.1
        self.u_R = 0.0
        
        # Adiabatic index
        self.gamma = 1.4
        
        # Discontinuity location
        self.x_diaphragm = 0.5
    
    def sound_speed(self, rho: float, p: float) -> float:
        """Calculate sound speed"""
        return np.sqrt(self.gamma * p / rho)
    
    def solve_riemann(self, t: float, x: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Solve Riemann problem analytically
        
        Args:
            t: Time
            x: Position array
            
        Returns:
            Dictionary with density, pressure, velocity, energy
        """
        # Sound speeds
        c_L = self.sound_speed(self.rho_L, self.p_L)
        c_R = self.sound_speed(self.rho_R, self.p_R)
        
        # Pressure ratio across shock
        p_ratio = self.p_L / self.p_R
        
        # Solve for post-shock pressure (iterative Newton-Raphson)
        p_star = self._find_p_star(p_ratio)
        
        # Post-shock velocity
        u_star = self._find_u_star(p_star)
        
        # Post-shock densities
        rho_star_L = self.rho_L * ((p_star / self.p_L + (self.gamma - 1) / (self.gamma + 1)) /
                                     ((self.gamma - 1) / (self.gamma + 1) * p_star / self.p_L + 1))
        rho_star_R = self.rho_R * (p_star / self.p_R) ** (1 / self.gamma)
        
        # Wave speeds
        shock_speed = u_star + c_R * np.sqrt((self.gamma + 1) / (2 * self.gamma) * p_star / self.p_R +
                                               (self.gamma - 1) / (2 * self.gamma))
        contact_speed = u_star
        rarefaction_head = -c_L
        rarefaction_tail = u_star - c_L * (p_star / self.p_L) ** ((self.gamma - 1) / (2 * self.gamma))
        
        # Initialize solution arrays
        rho = np.zeros_like(x)
        p = np.zeros_like(x)
        u = np.zeros_like(x)
        
        # Characteristic positions
        x_rel = (x - self.x_diaphragm) / t  # x/t characteristic
        
        for i, xi in enumerate(x_rel):
            if xi < rarefaction_head:
                # Left state
                rho[i] = self.rho_L
                p[i] = self.p_L
                u[i] = self.u_L
            elif xi < rarefaction_tail:
                # Rarefaction fan
                u[i] = 2 / (self.gamma + 1) * (c_L + xi)
                c = c_L - (self.gamma - 1) / 2 * u[i]
                rho[i] = self.rho_L * (c / c_L) ** (2 / (self.gamma - 1))
                p[i] = self.p_L * (c / c_L) ** (2 * self.gamma / (self.gamma - 1))
            elif xi < contact_speed:
                # Post-rarefaction state
                rho[i] = rho_star_L
                p[i] = p_star
                u[i] = u_star
            elif xi < shock_speed:
                # Post-shock state
                rho[i] = rho_star_R
                p[i] = p_star
                u[i] = u_star
            else:
                # Right state
                rho[i] = self.rho_R
                p[i] = self.p_R
                u[i] = self.u_R
        
        # Calculate energy (internal energy per unit mass)
        e = p / (rho * (self.gamma - 1))
        
        return {
            'density': rho,
            'pressure': p,
            'velocity': u,
            'energy': e
        }
    
    def _find_p_star(self, p_ratio: float, tol: float = 1e-6) -> float:
        """Find post-shock pressure using Newton-Raphson"""
        p_star = 0.5 * (self.p_L + self.p_R)  # Initial guess
        
        for _ in range(50):
            f = self._pressure_function(p_star)
            fp = self._pressure_derivative(p_star)
            p_new = p_star - f / fp
            
            if abs(p_new - p_star) < tol:
                return p_new
            p_star = p_new
        
        return p_star
    
    def _pressure_function(self, p: float) -> float:
        """Pressure function for Riemann solver"""
        c_L = self.sound_speed(self.rho_L, self.p_L)
        c_R = self.sound_speed(self.rho_R, self.p_R)
        
        # Left rarefaction
        f_L = 2 * c_L / (self.gamma - 1) * ((p / self.p_L) ** ((self.gamma - 1) / (2 * self.gamma)) - 1)
        
        # Right shock
        A_R = 2 / ((self.gamma + 1) * self.rho_R)
        B_R = (self.gamma - 1) / (self.gamma + 1) * self.p_R
        f_R = (p - self.p_R) * np.sqrt(A_R / (p + B_R))
        
        return f_L + f_R + self.u_R - self.u_L
    
    def _pressure_derivative(self, p: float) -> float:
        """Derivative of pressure function"""
        c_L = self.sound_speed(self.rho_L, self.p_L)
        
        df_L = (p / self.p_L) ** (-(self.gamma + 1) / (2 * self.gamma)) / (self.rho_L * c_L)
        
        A_R = 2 / ((self.gamma + 1) * self.rho_R)
        B_R = (self.gamma - 1) / (self.gamma + 1) * self.p_R
        df_R = np.sqrt(A_R / (p + B_R)) * (1 - (p - self.p_R) / (2 * (p + B_R)))
        
        return df_L + df_R
    
    def _find_u_star(self, p_star: float) -> float:
        """Calculate post-shock velocity"""
        c_L = self.sound_speed(self.rho_L, self.p_L)
        u_star = self.u_L - 2 * c_L / (self.gamma - 1) * ((p_star / self.p_L) ** ((self.gamma - 1) / (2 * self.gamma)) - 1)
        return u_star


def load_snapshot(filepath: Path) -> pd.DataFrame:
    """Load FDPS snapshot CSV file"""
    return pd.read_csv(filepath, comment='#')


def extract_1d_profile(data: pd.DataFrame, x_min: float = 0.0, x_max: float = 1.0, 
                        n_bins: int = 100) -> Dict[str, np.ndarray]:
    """
    Extract 1D profile from 3D SPH data by binning in x-direction
    
    Args:
        data: DataFrame with particle data
        x_min, x_max: Range in x-direction
        n_bins: Number of bins
        
    Returns:
        Dictionary with binned profiles
    """
    bins = np.linspace(x_min, x_max, n_bins + 1)
    x_centers = 0.5 * (bins[1:] + bins[:-1])
    
    # Bin particles
    indices = np.digitize(data['pos_x'], bins) - 1
    
    # Average properties in each bin
    density = np.zeros(n_bins)
    pressure = np.zeros(n_bins)
    velocity = np.zeros(n_bins)
    energy = np.zeros(n_bins)
    counts = np.zeros(n_bins)
    
    for i in range(n_bins):
        mask = indices == i
        if np.sum(mask) > 0:
            density[i] = data.loc[mask, 'dens'].mean()
            pressure[i] = data.loc[mask, 'pres'].mean()
            velocity[i] = data.loc[mask, 'vel_x'].mean()
            energy[i] = data.loc[mask, 'eng'].mean()
            counts[i] = np.sum(mask)
    
    # Handle empty bins by interpolation
    for arr in [density, pressure, velocity, energy]:
        mask = counts > 0
        if np.sum(mask) > 0:
            arr[counts == 0] = np.interp(x_centers[counts == 0], 
                                          x_centers[mask], 
                                          arr[mask])
    
    return {
        'x': x_centers,
        'density': density,
        'pressure': pressure,
        'velocity': velocity,
        'energy': energy,
        'counts': counts
    }


def plot_comparison(sim_profile: Dict, analytic_sol: Dict, time: float, 
                    output_file: Path):
    """
    Create comparison plot of simulation vs analytic solution
    
    Args:
        sim_profile: Simulation 1D profile
        analytic_sol: Analytic solution
        time: Simulation time
        output_file: Output file path
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'SPH Simulation vs Analytic Solution (t = {time:.4f})', 
                 fontsize=16, fontweight='bold')
    
    properties = [
        ('density', 'Density', 'g/cm³'),
        ('pressure', 'Pressure', 'dyne/cm²'),
        ('velocity', 'Velocity', 'cm/s'),
        ('energy', 'Specific Energy', 'erg/g')
    ]
    
    for ax, (prop, label, unit) in zip(axes.flat, properties):
        # Plot analytic solution
        ax.plot(analytic_sol['x'], analytic_sol[prop], 
                'k-', linewidth=2, label='Analytic', zorder=2)
        
        # Plot simulation
        ax.plot(sim_profile['x'], sim_profile[prop], 
                'ro', markersize=4, alpha=0.7, label='SPH Simulation', zorder=1)
        
        ax.set_xlabel('Position (x)', fontsize=11)
        ax.set_ylabel(f'{label} [{unit}]', fontsize=11)
        ax.legend(frameon=True, fancybox=True, shadow=True)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_title(label, fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✅ Created: {output_file}")
    plt.close()


def calculate_errors(sim_profile: Dict, analytic_sol: Dict) -> Dict[str, float]:
    """Calculate L1 and L2 errors"""
    errors = {}
    
    for prop in ['density', 'pressure', 'velocity', 'energy']:
        # Interpolate analytic solution to simulation grid
        analytic_interp = np.interp(sim_profile['x'], 
                                     analytic_sol['x'], 
                                     analytic_sol[prop])
        
        diff = sim_profile[prop] - analytic_interp
        
        # L1 error (mean absolute error)
        errors[f'{prop}_L1'] = np.mean(np.abs(diff))
        
        # L2 error (root mean square error)
        errors[f'{prop}_L2'] = np.sqrt(np.mean(diff ** 2))
        
        # Relative L2 error
        norm = np.sqrt(np.mean(analytic_interp ** 2))
        errors[f'{prop}_rel_L2'] = errors[f'{prop}_L2'] / (norm + 1e-10)
    
    return errors


def main():
    """Main comparison workflow"""
    print("=" * 70)
    print("FDPS SPH Simulation vs Analytic Solution Comparison")
    print("=" * 70)
    print()
    
    # Paths
    template_dir = Path(__file__).parent.parent
    result_dir = template_dir / "result"
    analysis_dir = template_dir / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    
    print(f"📁 Result directory: {result_dir}")
    print(f"📊 Analysis directory: {analysis_dir}")
    print()
    
    # Find all snapshots
    snapshots = sorted(result_dir.glob("snapshot_*_cgs.csv"))
    
    if not snapshots:
        print("❌ No snapshot files found!")
        print("   Please run the simulation first: make simulate")
        return 1
    
    print(f"📈 Found {len(snapshots)} snapshots")
    print()
    
    # Initialize Sod shock tube solver
    sod = SodShockTube()
    
    # Process each snapshot
    errors_history = []
    
    for i, snapshot_file in enumerate(snapshots):
        print(f"Processing {snapshot_file.name}...")
        
        # Load simulation data
        data = load_snapshot(snapshot_file)
        
        # Extract time from metadata
        with open(snapshot_file, 'r') as f:
            for line in f:
                if line.startswith('# Time:'):
                    time = float(line.split(':')[1].strip())
                    break
            else:
                time = i * 0.01  # Fallback
        
        # Extract 1D profile
        sim_profile = extract_1d_profile(data, x_min=0.0, x_max=1.0, n_bins=100)
        
        # Get analytic solution
        x_analytic = np.linspace(0.0, 1.0, 500)
        analytic_sol = sod.solve_riemann(time, x_analytic)
        analytic_sol['x'] = x_analytic
        
        # Create comparison plot
        output_file = analysis_dir / f"comparison_{i:04d}.png"
        plot_comparison(sim_profile, analytic_sol, time, output_file)
        
        # Calculate errors
        errors = calculate_errors(sim_profile, analytic_sol)
        errors['time'] = time
        errors['step'] = i
        errors_history.append(errors)
    
    # Save error history
    errors_df = pd.DataFrame(errors_history)
    errors_file = analysis_dir / "errors.csv"
    errors_df.to_csv(errors_file, index=False)
    print(f"\n✅ Error history saved: {errors_file}")
    
    # Plot error evolution
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Error Evolution Over Time', fontsize=16, fontweight='bold')
    
    properties = ['density', 'pressure', 'velocity', 'energy']
    labels = ['Density', 'Pressure', 'Velocity', 'Energy']
    
    for ax, prop, label in zip(axes.flat, properties, labels):
        ax.semilogy(errors_df['time'], errors_df[f'{prop}_L1'], 
                   'b-o', label='L1 Error', markersize=4)
        ax.semilogy(errors_df['time'], errors_df[f'{prop}_L2'], 
                   'r-s', label='L2 Error', markersize=4)
        ax.set_xlabel('Time', fontsize=11)
        ax.set_ylabel('Error', fontsize=11)
        ax.set_title(f'{label} Error', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    error_plot = analysis_dir / "error_evolution.png"
    plt.savefig(error_plot, dpi=150, bbox_inches='tight')
    print(f"✅ Created: {error_plot}")
    plt.close()
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("Summary Statistics (Final Snapshot)")
    print("=" * 70)
    final_errors = errors_history[-1]
    for prop in properties:
        print(f"\n{prop.capitalize()}:")
        print(f"  L1 Error:      {final_errors[f'{prop}_L1']:.6e}")
        print(f"  L2 Error:      {final_errors[f'{prop}_L2']:.6e}")
        print(f"  Relative L2:   {final_errors[f'{prop}_rel_L2']:.6e}")
    
    print("\n" + "=" * 70)
    print("✨ Comparison complete!")
    print(f"📂 All plots saved in: {analysis_dir}")
    print("=" * 70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
