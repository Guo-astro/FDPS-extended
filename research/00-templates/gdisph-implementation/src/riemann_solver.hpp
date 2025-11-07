#pragma once

/**
 * @file riemann_solver.hpp
 * @brief HLLC Riemann solver for GDISPH implementation
 * 
 * This implements the HLLC (Harten-Lax-van Leer-Contact) approximate Riemann solver
 * for compressible Euler equations. The solver computes numerical fluxes at particle
 * interfaces based on left and right primitive states.
 * 
 * References:
 * - GDISPH method: arXiv:2312.03224 (see /Users/guo-opt-p148/FDPS/papers/gdisph_2312.03224.pdf)
 * - HLLC solver: Toro, E.F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics.
 * 
 * TDD/BDD Implementation:
 * - RED phase: Created minimal structure ✓
 * - GREEN phase: Implemented to pass all tests ✓
 * - REFACTOR phase: Optimized with modern C++ practices ✓
 */

#include <particle_simulator.hpp>
#include <cmath>
#include <algorithm>

namespace GDISPH {

// Physical constants - use constexpr for compile-time evaluation
namespace Constants {
    constexpr PS::F64 kMinimumDensity = 1.0e-10;
    constexpr PS::F64 kMinimumPressure = 0.0;
    constexpr PS::F64 kSmallValue = 1.0e-10;  // For division guards
    constexpr PS::F64 kDefaultGamma = 1.4;     // Ideal diatomic gas
    constexpr PS::F64 kMinimumGamma = 1.0;     // Physical constraint
} // namespace Constants

/**
 * @brief Structure to hold computed flux values
 * 
 * For GDISPH (Godunov SPH), we solve the Riemann problem for internal energy only.
 * This is the correct formulation for Lagrangian methods where particles move with
 * the fluid. The kinetic energy is handled through particle advection, not flux computation.
 * 
 * Reference: arXiv:2312.03224 - GDISPH method
 */
struct RiemannFlux {
    PS::F64 mass_flux;          // Mass flux: rho * u_n
    PS::F64vec momentum_flux;   // Momentum flux: rho * v * u_n + P * n
    PS::F64 energy_flux;        // INTERNAL energy flux: (E_internal + P) * u_n
    PS::F64 pressure_star;      // Pressure in the star region (P*)
    
    RiemannFlux() : mass_flux(0.0), momentum_flux(0.0), energy_flux(0.0), pressure_star(0.0) {}
};

/**
 * @brief HLLC Riemann solver for compressible Euler equations
 * 
 * The HLLC solver is a three-wave solver that resolves contact discontinuities
 * better than the simpler HLL solver. It's particularly suited for SPH where
 * accurate interface resolution is crucial.
 */
class HLLCRiemannSolver {
private:
    PS::F64 gamma_;  // Adiabatic index (ratio of specific heats)
    
    /**
     * @brief Compute sound speed from primitive variables
     * @param rho Density (must be positive)
     * @param P Pressure (must be non-negative)
     * @return Sound speed, or 0 for invalid input
     * 
     * Invariant: Returns physically valid sound speed >= 0
     */
    constexpr PS::F64 computeSoundSpeed(PS::F64 rho, PS::F64 P) const noexcept {
        // c = sqrt(gamma * P / rho)
        // Guard against non-physical values using std::isfinite
        if (rho <= Constants::kMinimumDensity || P < Constants::kMinimumPressure || 
            !std::isfinite(rho) || !std::isfinite(P)) {
            return 0.0;
        }
        return std::sqrt(gamma_ * P / rho);
    }
    
    /**
     * @brief Compute specific internal energy from pressure
     * @param P Pressure
     * @param rho Density
     * @return Specific internal energy u = P/((gamma-1)*rho)
     * 
     * For Lagrangian SPH (GDISPH), we solve the Riemann problem using
     * internal energy only. Kinetic energy is handled separately through
     * particle advection in the Lagrangian frame.
     */
    constexpr PS::F64 computeInternalEnergy(PS::F64 P, PS::F64 rho) const noexcept {
        if (rho <= Constants::kMinimumDensity) {
            return 0.0;
        }
        return P / ((gamma_ - 1.0) * rho);
    }
    
    /**
     * @brief Estimate wave speeds using simple approach
     * 
     * For HLLC, we need three wave speeds:
     * - S_L: Left wave speed (fastest left-going wave)
     * - S_R: Right wave speed (fastest right-going wave)
     * - S_star: Contact wave speed (middle wave)
     */
    void estimateWaveSpeeds(PS::F64 rho_L, PS::F64 u_L, PS::F64 P_L, PS::F64 c_L,
                           PS::F64 rho_R, PS::F64 u_R, PS::F64 P_R, PS::F64 c_R,
                           PS::F64& S_L, PS::F64& S_R, PS::F64& S_star) const {
        // GREEN PHASE: Proper HLLC wave speed estimates
        // Using Davis direct wave speed estimates (simple and robust)
        S_L = std::min(u_L - c_L, u_R - c_R);
        S_R = std::max(u_L + c_L, u_R + c_R);
        
        // Contact wave speed (Equation 10.37 in Toro's book)
        PS::F64 rho_L_SL = rho_L * (S_L - u_L);
        PS::F64 rho_R_SR = rho_R * (S_R - u_R);
        
        PS::F64 numerator = P_R - P_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R);
        PS::F64 denominator = rho_L_SL - rho_R_SR;
        
        // Guard against division by zero
        if (std::abs(denominator) < 1.0e-10) {
            S_star = 0.5 * (u_L + u_R);
        } else {
            S_star = numerator / denominator;
        }
    }
    
    /**
     * @brief Compute conservative variables from primitive variables
     * 
     * For GDISPH (Lagrangian SPH with Godunov flux), we use:
     * - Mass density: rho
     * - Momentum density: rho * v
     * - Internal energy density: rho * u = P/(gamma-1)
     * 
     * Note: In Lagrangian SPH, kinetic energy is not part of the 
     * Riemann problem since particles move with the flow.
     */
    void computeConservativeVars(PS::F64 rho, const PS::F64vec& vel, PS::F64 P,
                                PS::F64& U_rho, PS::F64vec& U_mom, PS::F64& U_E_internal) const {
        U_rho = rho;
        U_mom = rho * vel;
        // Internal energy density (NOT total energy)
        U_E_internal = P / (gamma_ - 1.0);
    }
    
    /**
     * @brief Compute flux from conservative variables and velocity
     * 
     * For GDISPH, the flux computation uses internal energy:
     * - F_mass = rho * u_n
     * - F_momentum = rho * v * u_n + P * n
     * - F_internal_energy = (E_internal + P) * u_n
     * 
     * where E_internal = P/(gamma-1) is the internal energy density
     */
    void computeFluxFromConservative(PS::F64 U_rho, const PS::F64vec& U_mom, PS::F64 U_E_internal,
                                    PS::F64 u_n, PS::F64 P, const PS::F64vec& normal,
                                    PS::F64& F_rho, PS::F64vec& F_mom, PS::F64& F_E_internal) const {
        F_rho = U_rho * u_n;
        F_mom = U_mom * u_n;
        // Add pressure term in normal direction (momentum flux)
        F_mom = F_mom + P * normal;
        // Internal energy flux
        F_E_internal = (U_E_internal + P) * u_n;
    }
    
    /**
     * @brief Compute starred (intermediate) state for HLLC
     * 
     * This computes the star state in the HLLC solver for internal energy formulation.
     * The star region appears between the two wave speeds S_L and S_R.
     */
    void computeStarState(PS::F64 S_K, PS::F64 S_star,
                         PS::F64 U_rho_K, const PS::F64vec& U_mom_K, PS::F64 U_E_internal_K,
                         PS::F64 u_K, PS::F64 P_K,
                         PS::F64& U_rho_star, PS::F64vec& U_mom_star, PS::F64& U_E_internal_star) const {
        PS::F64 factor = U_rho_K * (S_K - u_K) / (S_K - S_star);
        
        // Density in star region
        U_rho_star = factor;
        
        // Momentum in star region (only normal component changes)
        U_mom_star = U_mom_K;
        U_mom_star.x = factor * S_star;  // Update normal component
        
        // Internal energy in star region
        // For internal energy formulation: E_star = E_K + (P_K / rho_K) * (S_star - u_K) / (S_K - u_K)
        PS::F64 u_internal_K = U_E_internal_K / U_rho_K;  // Specific internal energy
        U_E_internal_star = factor * (u_internal_K + (P_K / U_rho_K) * (S_star - u_K) / (S_K - u_K));
    }
    
public:
    /**
     * @brief Constructor
     * @param gamma Adiabatic index (default: 1.4 for ideal diatomic gas)
     */
    explicit HLLCRiemannSolver(PS::F64 gamma = 1.4) : gamma_(gamma) {
        // Validate gamma
        if (gamma <= 1.0) {
            std::cerr << "Warning: gamma must be > 1.0, using default 1.4" << std::endl;
            gamma_ = 1.4;
        }
    }
    
    /**
     * @brief Compute HLLC flux given left and right states
     * 
     * @param rho_L Left density
     * @param vel_L Left velocity vector
     * @param P_L Left pressure
     * @param rho_R Right density
     * @param vel_R Right velocity vector
     * @param P_R Right pressure
     * @param normal Unit normal vector pointing from left to right
     * @return RiemannFlux structure containing mass, momentum, and energy fluxes
     * 
     * GIVEN: Left and right primitive states (rho, vel, P) and interface normal
     * WHEN: computeFlux is called
     * THEN: Returns physically consistent flux values
     */
    RiemannFlux computeFlux(PS::F64 rho_L, const PS::F64vec& vel_L, PS::F64 P_L,
                           PS::F64 rho_R, const PS::F64vec& vel_R, PS::F64 P_R,
                           const PS::F64vec& normal) const {
        RiemannFlux flux;
        
        // GREEN PHASE: Full HLLC implementation
        
        // Step 1: Rotate velocities to align with normal (1D Riemann problem)
        PS::F64 u_L = vel_L * normal;  // Normal velocity component
        PS::F64 u_R = vel_R * normal;
        
        // Step 2: Compute sound speeds
        PS::F64 c_L = computeSoundSpeed(rho_L, P_L);
        PS::F64 c_R = computeSoundSpeed(rho_R, P_R);
        
        // Step 3: Estimate wave speeds
        PS::F64 S_L, S_R, S_star;
        estimateWaveSpeeds(rho_L, u_L, P_L, c_L, rho_R, u_R, P_R, c_R, S_L, S_R, S_star);
        
        // Step 4: Compute P_star (pressure in the star region)
        // From HLLC: P_star = P_K + rho_K * (S_K - u_K) * (S_star - u_K)
        // Should be the same from left or right; we average for robustness
        PS::F64 P_star_L = P_L + rho_L * (S_L - u_L) * (S_star - u_L);
        PS::F64 P_star_R = P_R + rho_R * (S_R - u_R) * (S_star - u_R);
        PS::F64 P_star = 0.5 * (P_star_L + P_star_R);
        
        // Store P_star in flux structure for GDISPH
        flux.pressure_star = P_star;
        
        // Step 5: Compute fluxes based on wave configuration
        
        if (S_L >= 0.0) {
            // Supersonic flow from left
            PS::F64 U_rho_L, U_E_internal_L;
            PS::F64vec U_mom_L;
            computeConservativeVars(rho_L, vel_L, P_L, U_rho_L, U_mom_L, U_E_internal_L);
            
            PS::F64 F_rho, F_E_internal;
            PS::F64vec F_mom;
            computeFluxFromConservative(U_rho_L, U_mom_L, U_E_internal_L, u_L, P_L, normal, 
                                       F_rho, F_mom, F_E_internal);
            
            flux.mass_flux = F_rho;
            flux.momentum_flux = F_mom;
            flux.energy_flux = F_E_internal;
            
        } else if (S_R <= 0.0) {
            // Supersonic flow from right
            PS::F64 U_rho_R, U_E_internal_R;
            PS::F64vec U_mom_R;
            computeConservativeVars(rho_R, vel_R, P_R, U_rho_R, U_mom_R, U_E_internal_R);
            
            PS::F64 F_rho, F_E_internal;
            PS::F64vec F_mom;
            computeFluxFromConservative(U_rho_R, U_mom_R, U_E_internal_R, u_R, P_R, normal, 
                                       F_rho, F_mom, F_E_internal);
            
            flux.mass_flux = F_rho;
            flux.momentum_flux = F_mom;
            flux.energy_flux = F_E_internal;
            
        } else if (S_L < 0.0 && S_star >= 0.0) {
            // Left star state
            PS::F64 U_rho_L, U_E_internal_L;
            PS::F64vec U_mom_L;
            computeConservativeVars(rho_L, vel_L, P_L, U_rho_L, U_mom_L, U_E_internal_L);
            
            // Compute left flux
            PS::F64 F_rho_L, F_E_internal_L;
            PS::F64vec F_mom_L;
            computeFluxFromConservative(U_rho_L, U_mom_L, U_E_internal_L, u_L, P_L, normal, 
                                       F_rho_L, F_mom_L, F_E_internal_L);
            
            // Compute star state
            PS::F64 U_rho_star, U_E_internal_star;
            PS::F64vec U_mom_star;
            computeStarState(S_L, S_star, U_rho_L, U_mom_L, U_E_internal_L, u_L, P_L,
                           U_rho_star, U_mom_star, U_E_internal_star);
            
            // HLLC flux formula: F_star = F_L + S_L * (U_star - U_L)
            flux.mass_flux = F_rho_L + S_L * (U_rho_star - U_rho_L);
            flux.momentum_flux = F_mom_L + S_L * (U_mom_star - U_mom_L);
            flux.energy_flux = F_E_internal_L + S_L * (U_E_internal_star - U_E_internal_L);
            
        } else {
            // Right star state (S_star < 0.0 && S_R > 0.0)
            PS::F64 U_rho_R, U_E_internal_R;
            PS::F64vec U_mom_R;
            computeConservativeVars(rho_R, vel_R, P_R, U_rho_R, U_mom_R, U_E_internal_R);
            
            // Compute right flux
            PS::F64 F_rho_R, F_E_internal_R;
            PS::F64vec F_mom_R;
            computeFluxFromConservative(U_rho_R, U_mom_R, U_E_internal_R, u_R, P_R, normal, 
                                       F_rho_R, F_mom_R, F_E_internal_R);
            
            // Compute star state
            PS::F64 U_rho_star, U_E_internal_star;
            PS::F64vec U_mom_star;
            computeStarState(S_R, S_star, U_rho_R, U_mom_R, U_E_internal_R, u_R, P_R,
                           U_rho_star, U_mom_star, U_E_internal_star);
            
            // HLLC flux formula: F_star = F_R + S_R * (U_star - U_R)
            flux.mass_flux = F_rho_R + S_R * (U_rho_star - U_rho_R);
            flux.momentum_flux = F_mom_R + S_R * (U_mom_star - U_mom_R);
            flux.energy_flux = F_E_internal_R + S_R * (U_E_internal_star - U_E_internal_R);
        }
        
        return flux;
    }
    
    /**
     * @brief Get the adiabatic index
     */
    PS::F64 getGamma() const { return gamma_; }
};

} // namespace GDISPH
