#pragma once

#include <particle_simulator.hpp>
#include <cmath>
#include <iostream>
#include <cassert>
#include <string>

// Include the Riemann solver we'll implement
#include "../src/riemann_solver.hpp"

namespace GDISPH {
namespace Test {

/**
 * @brief Test helper to compare floating point values with tolerance
 */
constexpr PS::F64 kTestTolerance = 1.0e-10;

inline bool isClose(PS::F64 a, PS::F64 b, PS::F64 tol = kTestTolerance) {
    return std::abs(a - b) < tol;
}

inline bool isClose(const PS::F64vec& a, const PS::F64vec& b, PS::F64 tol = kTestTolerance) {
    return std::abs(a.x - b.x) < tol && 
           std::abs(a.y - b.y) < tol && 
           std::abs(a.z - b.z) < tol;
}

/**
 * @brief Test result tracker
 */
struct TestResult {
    std::string test_name;
    bool passed{true};
    std::string failure_message;
};

class RiemannSolverTestSuite {
private:
    std::vector<TestResult> results_;
    
    void recordTest(const std::string& name, bool passed, const std::string& msg = "") {
        results_.push_back({name, passed, msg});
        if (!passed) {
            std::cerr << "❌ FAILED: " << name << std::endl;
            if (!msg.empty()) {
                std::cerr << "   " << msg << std::endl;
            }
        } else {
            std::cout << "✅ PASSED: " << name << std::endl;
        }
    }
    
public:
    /**
     * GIVEN: Two identical states at rest
     * WHEN: HLLC solver computes the interface flux
     * THEN: Flux should have zero mass and energy flux, pressure in momentum
     */
    void test_GivenIdenticalStatesAtRest_WhenComputingFlux_ThenFluxIsZero() {
        // Given: Identical states at rest
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 1.0;
        PS::F64vec vel_L(0.0, 0.0, 0.0);
        PS::F64vec vel_R(0.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 1.0;
        PS::F64 gamma = 1.4;
        PS::F64vec normal(1.0, 0.0, 0.0);
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: Mass flux and energy flux should be zero (no flow)
        // Momentum flux contains pressure term which is non-zero but symmetric
        bool passed = isClose(flux.mass_flux, 0.0, 1e-10) &&
                     isClose(flux.energy_flux, 0.0, 1e-10);
        
        recordTest("GivenIdenticalStatesAtRest_WhenComputingFlux_ThenFluxIsZero", passed,
                  passed ? "" : "Expected zero mass and energy flux for identical states at rest");
    }
    
    /**
     * GIVEN: Left state moving with constant velocity, right state at rest
     * WHEN: HLLC solver computes the interface flux
     * THEN: Flux should capture mass and momentum transport
     */
    void test_GivenLeftStateMoving_WhenComputingFlux_ThenCapturesTransport() {
        // Given: Left moving, right at rest
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 1.0;
        PS::F64vec vel_L(1.0, 0.0, 0.0);
        PS::F64vec vel_R(0.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 1.0;
        PS::F64 gamma = 1.4;
        PS::F64vec normal(1.0, 0.0, 0.0);
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: Should have non-zero mass flux
        bool passed = flux.mass_flux > 0.0;
        
        recordTest("GivenLeftStateMoving_WhenComputingFlux_ThenCapturesTransport", passed,
                  passed ? "" : "Expected positive mass flux for left state moving right");
    }
    
    /**
     * GIVEN: Sod shock tube initial conditions (classic test)
     * WHEN: HLLC solver computes the interface flux
     * THEN: Flux should be physically reasonable (positive wave speeds)
     */
    void test_GivenSodShockTubeIC_WhenComputingFlux_ThenFluxIsPhysical() {
        // Given: Sod shock tube (high pressure left, low pressure right)
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 0.125;
        PS::F64vec vel_L(0.0, 0.0, 0.0);
        PS::F64vec vel_R(0.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 0.1;
        PS::F64 gamma = 1.4;
        PS::F64vec normal(1.0, 0.0, 0.0);
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: Flux should be non-zero and finite
        bool passed = std::isfinite(flux.mass_flux) && 
                     std::isfinite(flux.momentum_flux.x) &&
                     std::isfinite(flux.energy_flux);
        
        recordTest("GivenSodShockTubeIC_WhenComputingFlux_ThenFluxIsPhysical", passed,
                  passed ? "" : "Flux values must be finite for Sod shock tube");
    }
    
    /**
     * GIVEN: States with different normal directions
     * WHEN: HLLC solver computes flux
     * THEN: Flux should be correctly projected along normal
     */
    void test_GivenDifferentNormals_WhenComputingFlux_ThenFluxProjectedCorrectly() {
        // Given: States with y-direction normal
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 1.0;
        PS::F64vec vel_L(0.0, 1.0, 0.0);  // Moving in y
        PS::F64vec vel_R(0.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 1.0;
        PS::F64 gamma = 1.4;
        PS::F64vec normal(0.0, 1.0, 0.0);  // Normal in y direction
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: Should have positive mass flux (flow along normal)
        bool passed = flux.mass_flux > 0.0;
        
        recordTest("GivenDifferentNormals_WhenComputingFlux_ThenFluxProjectedCorrectly", passed,
                  passed ? "" : "Expected positive mass flux along y-normal");
    }
    
    /**
     * GIVEN: High Mach number flow (supersonic)
     * WHEN: HLLC solver computes flux
     * THEN: Solver should remain stable and return finite values
     */
    void test_GivenSupersonicFlow_WhenComputingFlux_ThenSolverIsStable() {
        // Given: High velocity (Mach >> 1)
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 0.1;
        PS::F64vec vel_L(10.0, 0.0, 0.0);  // High velocity
        PS::F64vec vel_R(0.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 0.1;
        PS::F64 gamma = 1.4;
        PS::F64vec normal(1.0, 0.0, 0.0);
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: All values should be finite (no NaN/Inf)
        bool passed = std::isfinite(flux.mass_flux) && 
                     std::isfinite(flux.momentum_flux.x) &&
                     std::isfinite(flux.momentum_flux.y) &&
                     std::isfinite(flux.momentum_flux.z) &&
                     std::isfinite(flux.energy_flux);
        
        recordTest("GivenSupersonicFlow_WhenComputingFlux_ThenSolverIsStable", passed,
                  passed ? "" : "Solver must remain stable for supersonic flows");
    }
    
    /**
     * GIVEN: Vacuum-like state (very low density/pressure)
     * WHEN: HLLC solver computes flux
     * THEN: Solver should handle near-vacuum without crashing
     */
    void test_GivenNearVacuum_WhenComputingFlux_ThenHandlesGracefully() {
        // Given: Near-vacuum state
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 1.0e-10;  // Very low density
        PS::F64vec vel_L(0.0, 0.0, 0.0);
        PS::F64vec vel_R(0.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 1.0e-10;  // Very low pressure
        PS::F64 gamma = 1.4;
        PS::F64vec normal(1.0, 0.0, 0.0);
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: Should not crash and return finite values
        bool passed = std::isfinite(flux.mass_flux) && 
                     std::isfinite(flux.momentum_flux.x) &&
                     std::isfinite(flux.energy_flux);
        
        recordTest("GivenNearVacuum_WhenComputingFlux_ThenHandlesGracefully", passed,
                  passed ? "" : "Solver must handle near-vacuum states gracefully");
    }
    
    /**
     * GIVEN: Symmetric shock problem (same magnitude, opposite velocities)
     * WHEN: HLLC solver computes flux
     * THEN: Flux should show symmetry properties
     */
    void test_GivenSymmetricShock_WhenComputingFlux_ThenShowsSymmetry() {
        // Given: Symmetric collision
        PS::F64 rho_L = 1.0;
        PS::F64 rho_R = 1.0;
        PS::F64vec vel_L(1.0, 0.0, 0.0);
        PS::F64vec vel_R(-1.0, 0.0, 0.0);
        PS::F64 P_L = 1.0;
        PS::F64 P_R = 1.0;
        PS::F64 gamma = 1.4;
        PS::F64vec normal(1.0, 0.0, 0.0);
        
        // When: Computing flux
        HLLCRiemannSolver solver(gamma);
        auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
        
        // Then: Mass flux should be zero (symmetric)
        bool passed = isClose(flux.mass_flux, 0.0, 1e-8);
        
        recordTest("GivenSymmetricShock_WhenComputingFlux_ThenShowsSymmetry", passed,
                  passed ? "" : "Expected zero mass flux for symmetric collision");
    }
    
    /**
     * Run all tests and report summary
     */
    void runAllTests() {
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "Running Riemann Solver Test Suite (BDD Style)" << std::endl;
        std::cout << std::string(60, '=') << "\n" << std::endl;
        
        test_GivenIdenticalStatesAtRest_WhenComputingFlux_ThenFluxIsZero();
        test_GivenLeftStateMoving_WhenComputingFlux_ThenCapturesTransport();
        test_GivenSodShockTubeIC_WhenComputingFlux_ThenFluxIsPhysical();
        test_GivenDifferentNormals_WhenComputingFlux_ThenFluxProjectedCorrectly();
        test_GivenSupersonicFlow_WhenComputingFlux_ThenSolverIsStable();
        test_GivenNearVacuum_WhenComputingFlux_ThenHandlesGracefully();
        test_GivenSymmetricShock_WhenComputingFlux_ThenShowsSymmetry();
        
        // Summary
        int passed = 0;
        int failed = 0;
        for (const auto& result : results_) {
            if (result.passed) passed++;
            else failed++;
        }
        
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "Test Summary:" << std::endl;
        std::cout << "  ✅ Passed: " << passed << std::endl;
        std::cout << "  ❌ Failed: " << failed << std::endl;
        std::cout << "  Total:  " << (passed + failed) << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        
        if (failed > 0) {
            std::cout << "\n⚠️  Some tests failed! Review failures above." << std::endl;
        } else {
            std::cout << "\n🎉 All tests passed!" << std::endl;
        }
    }
    
    bool allTestsPassed() const {
        for (const auto& result : results_) {
            if (!result.passed) return false;
        }
        return true;
    }
};

} // namespace Test
} // namespace GDISPH
