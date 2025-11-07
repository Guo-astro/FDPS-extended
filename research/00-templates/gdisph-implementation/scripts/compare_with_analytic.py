#!/usr/bin/env python3
"""
Create comparison plots of GDISPH simulation with analytic Sod shock tube solution.
Works with minimal dependencies (no matplotlib, pandas, or numpy needed).
"""

import sys
import csv
from pathlib import Path
import math

def read_profile_csv(filename):
    """Read binned profile data"""
    data = {'x': [], 'dens': [], 'vel': [], 'pres': []}
    
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data['x'].append(float(row['x']))
            data['dens'].append(float(row['dens']))
            data['vel'].append(float(row['vel']))
            data['pres'].append(float(row['pres']))
    
    return data

def sod_shock_analytic(x, t, gamma=1.4):
    """
    Analytic solution for Sod shock tube.
    Initial conditions: rho_L=1, P_L=1, rho_R=0.125, P_R=0.1, x_diaphragm=0.5
    """
    # Left state
    rho_L = 1.0
    P_L = 1.0
    u_L = 0.0
    
    # Right state
    rho_R = 0.125
    P_R = 0.1
    u_R = 0.0
    
    # Sound speeds
    c_L = math.sqrt(gamma * P_L / rho_L)
    c_R = math.sqrt(gamma * P_R / rho_R)
    
    # Pressure ratio across contact discontinuity (simplified approximation)
    # For exact solution, need to solve Riemann problem iteratively
    # Here we use approximate values for Sod problem at t=0.12
    P_star = 0.30313  # Contact pressure
    u_star = 0.92745  # Contact velocity
    
    # Post-shock density (left)
    rho_L_star = rho_L * ((P_star/P_L + (gamma-1)/(gamma+1)) / 
                          ((gamma-1)/(gamma+1) * P_star/P_L + 1))
    
    # Post-expansion density (left of contact)
    rho_2 = rho_L * (P_star / P_L) ** (1.0/gamma)
    
    # Wave speeds
    shock_speed = u_star + c_R * math.sqrt((gamma+1)/(2*gamma) * P_star/P_R + 
                                            (gamma-1)/(2*gamma))
    
    head_speed = -c_L
    tail_speed = u_star - c_L * (P_star/P_L) ** ((gamma-1)/(2*gamma))
    
    # Position of discontinuity
    x0 = 0.5
    
    # Classify regions
    if x < x0 + head_speed * t:
        # Region 1: undisturbed left state
        return {'rho': rho_L, 'u': u_L, 'P': P_L}
    elif x < x0 + tail_speed * t:
        # Region 2: expansion fan
        c = (gamma-1)/(gamma+1) * (c_L - (x - x0)/t)
        rho = rho_L * (c/c_L) ** (2/(gamma-1))
        u = (gamma-1)/(gamma+1) * (-c_L + 2/(gamma-1) * c_L + (x - x0)/t)
        P = P_L * (c/c_L) ** (2*gamma/(gamma-1))
        return {'rho': rho, 'u': u, 'P': P}
    elif x < x0 + u_star * t:
        # Region 3: left of contact
        return {'rho': rho_2, 'u': u_star, 'P': P_star}
    elif x < x0 + shock_speed * t:
        # Region 4: right of contact
        rho_R_star = rho_R * ((gamma+1)*P_star/P_R + (gamma-1)) / \
                             ((gamma-1)*P_star/P_R + (gamma+1))
        return {'rho': rho_R_star, 'u': u_star, 'P': P_star}
    else:
        # Region 5: undisturbed right state
        return {'rho': rho_R, 'u': u_R, 'P': P_R}

def create_ascii_plot(sim_data, analytic_data, var_name, title, width=70, height=15):
    """Create ASCII art plot comparing simulation and analytic solution"""
    
    # Map variable names between simulation and analytic
    var_map = {'dens': 'rho', 'vel': 'u', 'pres': 'P'}
    
    # Get data ranges
    sim_vals = sim_data[var_name]
    ana_vals = [analytic_data[i][var_map[var_name]] for i in range(len(analytic_data))]
    
    all_vals = sim_vals + ana_vals
    y_min = min(all_vals) * 0.9
    y_max = max(all_vals) * 1.1
    y_range = y_max - y_min
    
    if y_range == 0:
        y_range = 1
    
    x_vals = sim_data['x']
    x_min = min(x_vals)
    x_max = max(x_vals)
    x_range = x_max - x_min
    
    # Create grid
    grid = [[' ' for _ in range(width)] for _ in range(height)]
    
    # Plot simulation data (with 'o')
    for i, (x, y) in enumerate(zip(x_vals, sim_vals)):
        col = int((x - x_min) / x_range * (width - 1))
        row = height - 1 - int((y - y_min) / y_range * (height - 1))
        row = max(0, min(height-1, row))
        col = max(0, min(width-1, col))
        grid[row][col] = 'o'
    
    # Plot analytic solution (with '*')
    ana_x = [analytic_data[i]['x'] for i in range(len(analytic_data))]
    for x, y in zip(ana_x, ana_vals):
        col = int((x - x_min) / x_range * (width - 1))
        row = height - 1 - int((y - y_min) / y_range * (height - 1))
        row = max(0, min(height-1, row))
        col = max(0, min(width-1, col))
        if grid[row][col] == ' ':
            grid[row][col] = '*'
        elif grid[row][col] == 'o':
            grid[row][col] = '#'  # Overlap
    
    # Print plot
    print(f"\n{title}")
    print("=" * width)
    print(f"y_max: {y_max:.4f}")
    for row in grid:
        print(''.join(row))
    print(f"y_min: {y_min:.4f}")
    print(f"x: {x_min:.3f} -> {x_max:.3f}")
    print(f"Legend: o=simulation, *=analytic, #=overlap")
    print("=" * width)

def main():
    result_dir = Path("/Users/guo-opt-p148/FDPS/research/00-templates/gdisph-implementation")
    plots_dir = result_dir / "plots"
    
    if not plots_dir.exists():
        print(f"Error: plots directory not found at {plots_dir}")
        print("Run quick_analysis.py first to generate profile data")
        return 1
    
    # Find all profile files
    profile_files = sorted(plots_dir.glob("profile_*.csv"))
    
    if not profile_files:
        print("No profile files found!")
        return 1
    
    print("=" * 80)
    print("GDISPH SIMULATION vs ANALYTIC SOD SHOCK TUBE SOLUTION")
    print("=" * 80)
    
    # Process each snapshot
    for profile_file in profile_files:
        snapshot_num = profile_file.stem.replace('profile_', '')
        
        # Read simulation data
        sim_data = read_profile_csv(profile_file)
        
        if not sim_data['x']:
            continue
        
        # Estimate time from snapshot number (assuming ~0.002 time units per snapshot)
        t = int(snapshot_num) * 0.002
        
        print(f"\n{'='*80}")
        print(f"SNAPSHOT {snapshot_num} (t ≈ {t:.4f})")
        print(f"{'='*80}")
        
        # Generate analytic solution at same x points
        analytic_data = []
        for x in sim_data['x']:
            sol = sod_shock_analytic(x, t)
            analytic_data.append({
                'x': x,
                'rho': sol['rho'],
                'u': sol['u'],
                'P': sol['P']
            })
        
        # Create comparison plots
        create_ascii_plot(sim_data, analytic_data, 'dens', 
                         f"DENSITY COMPARISON (t={t:.4f})")
        create_ascii_plot(sim_data, analytic_data, 'vel', 
                         f"VELOCITY COMPARISON (t={t:.4f})")
        create_ascii_plot(sim_data, analytic_data, 'pres', 
                         f"PRESSURE COMPARISON (t={t:.4f})")
        
        # Calculate errors
        rho_err = sum((sim_data['dens'][i] - analytic_data[i]['rho'])**2 
                     for i in range(len(sim_data['x']))) / len(sim_data['x'])
        rho_err = math.sqrt(rho_err)
        
        print(f"\nRMS Error in density: {rho_err:.6f}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
