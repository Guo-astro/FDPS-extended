#include <particle_simulator.hpp>
#include "../src/riemann_solver.hpp"
#include <iostream>
#include <cmath>

int main() {
    PS::F64 rho_L = 1.0;
    PS::F64 rho_R = 1.0;
    PS::F64vec vel_L(0.0, 0.0, 0.0);
    PS::F64vec vel_R(0.0, 0.0, 0.0);
    PS::F64 P_L = 1.0;
    PS::F64 P_R = 1.0;
    PS::F64 gamma = 1.4;
    PS::F64vec normal(1.0, 0.0, 0.0);
    
    GDISPH::HLLCRiemannSolver solver(gamma);
    auto flux = solver.computeFlux(rho_L, vel_L, P_L, rho_R, vel_R, P_R, normal);
    
    std::cout << "Identical states at rest test:" << std::endl;
    std::cout << "  mass_flux = " << flux.mass_flux << std::endl;
    std::cout << "  momentum_flux = (" << flux.momentum_flux.x << ", " 
              << flux.momentum_flux.y << ", " << flux.momentum_flux.z << ")" << std::endl;
    std::cout << "  energy_flux = " << flux.energy_flux << std::endl;
    
    return 0;
}
