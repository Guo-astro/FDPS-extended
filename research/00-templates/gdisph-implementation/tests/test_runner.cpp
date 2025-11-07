#include <iostream>
#include "test_riemann_solver.hpp"

int main(int argc, char* argv[]) {
    std::cout << "GDISPH Riemann Solver - Test-Driven Development" << std::endl;
    std::cout << "Following BDD (Behavior-Driven Development) approach" << std::endl;
    std::cout << std::endl;
    
    GDISPH::Test::RiemannSolverTestSuite test_suite;
    test_suite.runAllTests();
    
    return test_suite.allTestsPassed() ? 0 : 1;
}
