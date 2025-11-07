# GDISPH Implementation Notes

## Overview

This directory contains an implementation of **GDISPH (Godunov Smoothed Particle Hydrodynamics)** using the FDPS framework. GDISPH is an advanced SPH method that uses Riemann solvers for improved shock capturing and conservation properties.

## Key References

- **Primary Reference**: arXiv:2312.03224 - "Novel Hydrodynamic Schemes Capturing Shocks and Contact Discontinuities and Comparison Study with Existing Methods"
  - Paper location: `/Users/guo-opt-p148/FDPS/papers/gdisph_2312.03224.pdf`
- **HLLC Riemann Solver**: Toro, E.F. (2009). "Riemann Solvers and Numerical Methods for Fluid Dynamics"

**Note**: While Inutsuka (2002) introduced the concept of using Riemann solvers in SPH (GSPH), the specific GDISPH method implemented here follows the formulation in arXiv:2312.03224.

## Implementation Details

### 1. Riemann Solver (`src/riemann_solver.hpp`)

Implements the **HLLC (Harten-Lax-van Leer-Contact)** Riemann solver for compressible Euler equations.

**Key Design Decision: Internal Energy Only**

For Lagrangian SPH methods like GDISPH:
- ✅ **Use internal energy** in Riemann solver
- ❌ **NOT total energy** (internal + kinetic)

**Rationale:**
- SPH particles move with the fluid (Lagrangian frame)
- Kinetic energy is handled through particle advection
- Only internal (thermal) energy requires flux computation
- This is the correct formulation for all Lagrangian hydro methods

**Implementation:**
```cpp
// Conservative variables
U_rho = rho
U_mom = rho * v  
U_E_internal = P/(gamma-1)  // NOT P/(gamma-1) + 0.5*rho*v^2

// Fluxes
F_mass = rho * u_n
F_momentum = rho * v * u_n + P * n
F_energy = (E_internal + P) * u_n  // Internal energy flux only
```

**Test Coverage:**
- 7/7 BDD-style tests passing
- Tests cover: identical states, moving states, Sod shock tube, different normals, supersonic flow, near-vacuum, symmetric shocks

### 2. GDISPH Force Calculation (`src/gdisph_force.hpp`)

Implements force and energy updates using Riemann fluxes.

**Momentum Equation** (arXiv:2312.03224):
```
dv_i/dt = -sum_j m_j * (P*_ij/rho_i^2 + P*_ji/rho_j^2) * grad_i W_ij
```

Where:
- `P*_ij` = contact pressure from Riemann solver between particles i and j
- `grad_i W_ij` = kernel gradient

**Energy Equation** (arXiv:2312.03224):
```
du_i/dt = sum_j m_j * (P*_ij/rho_i^2) * (v*_ij - v_i) · grad_i W_ij
```

Where:
- `u` = specific internal energy
- `v*_ij` = contact velocity from Riemann solver

**No Artificial Viscosity:**
- GDISPH does NOT use artificial viscosity
- Shock capturing is handled naturally by the Riemann solver
- This is a major advantage over standard SPH

### 3. Timestep Calculation

Signal velocity for CFL condition:
```
v_sig = v_approach + c_i + c_j
dt = CFL * 2*h / v_sig_max
```

Where:
- `v_approach` = max(0, -(v_i - v_j) · n_ij) (approaching velocity)
- `c_i`, `c_j` = sound speeds
- `CFL` = Courant number (typically 0.3)

## Differences from Standard SPH

| Feature | Standard SPH | GDISPH |
|---------|-------------|--------|
| **Pressure force** | Direct pressure gradient | Riemann flux (contact pressure P*) |
| **Shock handling** | Artificial viscosity | Natural (via Riemann solver) |
| **Energy** | Total energy or internal | Internal energy only (Lagrangian) |
| **Conservation** | Approximate | Improved (Godunov flux) |
| **Contact discontinuities** | Smeared | Better resolved (HLLC) |

## File Structure

```
gdisph-implementation/
├── src/
│   ├── main.cpp                 # Main simulation (uses GDISPH)
│   ├── riemann_solver.hpp       # HLLC Riemann solver
│   └── gdisph_force.hpp         # GDISPH force calculation
├── tests/
│   ├── test_riemann_solver.hpp  # BDD-style test suite
│   ├── test_runner.cpp          # Test executable
│   └── Makefile                 # Test build system
└── GDISPH_IMPLEMENTATION_NOTES.md  # This file
```

## Building and Testing

### Test the Riemann Solver

```bash
cd tests
nix develop /path/to/FDPS --command bash -c 'make test'
```

Expected output: **7/7 tests passing** ✅

### Build Main Simulation

```bash
nix develop /path/to/FDPS --command bash -c 'make all'
```

### Run Simulation

```bash
./scripts/run_simulation.sh
```

## TDD/BDD Development Process

This implementation followed Test-Driven Development with Behavior-Driven Design:

1. **RED Phase**: Created failing test suite
2. **GREEN Phase**: Implemented HLLC solver to pass tests
3. **REFACTOR Phase**: Optimized for modern C++ (constexpr, noexcept)

All test names follow Given-When-Then pattern:
```cpp
test_GivenIdenticalStatesAtRest_WhenComputingFlux_ThenFluxIsZero()
test_GivenSodShockTubeIC_WhenComputingFlux_ThenFluxIsPhysical()
// etc.
```

## Next Steps / TODO

- [ ] Add second-order reconstruction (gradient estimation + slope limiting)
- [ ] Implement Sod shock tube validation test
- [ ] Compare with standard SPH on shock problems
- [ ] Add more test cases (blast wave, Sedov, etc.)
- [ ] Performance profiling and optimization
- [ ] Document expected results and convergence rates

## Common Issues / FAQ

**Q: Why does momentum flux include pressure term?**

A: The momentum flux is F_mom = ρv·u_n + P·n. The pressure term gives the force contribution.

**Q: Do we still need artificial viscosity?**

A: No! The Riemann solver handles shocks naturally. This is a key advantage of GDISPH.

**Q: What about tensile instability?**

A: GDISPH still inherits SPH's tensile instability in tension. May need additional stabilization for negative pressure states.

**Q: Can I use total energy instead?**

A: No! For Lagrangian methods, you must use internal energy only. Total energy would double-count kinetic energy.

## Validation Recommendations

1. **Sod Shock Tube** (1D): Standard test for shock-capturing methods
2. **Sedov Blast Wave** (2D/3D): Strong shock propagation
3. **Kelvin-Helmholtz Instability**: Test for contact discontinuities
4. **Interacting Blast Waves**: Multiple shock interactions

Compare results with:
- Exact Riemann solutions (where available)
- High-resolution Godunov finite-volume codes
- Standard SPH with artificial viscosity

## Contact & Support

For questions about this implementation, refer to:
- FDPS documentation: /doc/
- GDISPH paper: arXiv:2312.03224 (local: /Users/guo-opt-p148/FDPS/papers/gdisph_2312.03224.pdf)
- HLLC solver: Toro (2009) "Riemann Solvers and Numerical Methods for Fluid Dynamics"
