/*
 * Unit tests for FDPS Unit System and Output Formats
 * 
 * This file demonstrates:
 * 1. Compile-time dimension checking
 * 2. Unit conversions between different systems
 * 3. Output format functionality
 * 4. Extensibility of the system
 */

#include <iostream>
#include <cassert>
#include <cmath>
#include <particle_simulator.hpp>
#include "fdps_unit_system.hpp"
#include "fdps_output.hpp"

using namespace PS;
using namespace PS::UnitSystem;
using namespace PS::Output;

// Test particle structure
struct TestParticle {
    S64 id;
    F64 mass;
    F64vec pos;
    F64vec vel;
    F64 dens;
    F64 eng;
    F64 pres;
    F64 smth;
    F64 snds;
    F64 eng_dot;
    F64 dt;
    F64vec acc;
    F64vec vel_half;
    F64 eng_half;
};

// Helper function for floating point comparison
bool isClose(F64 a, F64 b, F64 rel_tol = 1e-9, F64 abs_tol = 1e-12) {
    return std::abs(a - b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol);
}

// ===========================================
// Test 1: Basic quantity creation and arithmetic
// ===========================================
void test_quantity_basic() {
    std::cout << "Test 1: Basic quantity operations..." << std::endl;
    
    Mass<CGSUnit> m1(100.0);
    Mass<CGSUnit> m2(50.0);
    
    // Addition
    Mass<CGSUnit> m_sum = m1 + m2;
    assert(m_sum.getValue() == 150.0);
    
    // Subtraction
    Mass<CGSUnit> m_diff = m1 - m2;
    assert(m_diff.getValue() == 50.0);
    
    // Scalar multiplication
    Mass<CGSUnit> m_scaled = m1 * 2.0;
    assert(m_scaled.getValue() == 200.0);
    
    // Scalar division
    Mass<CGSUnit> m_divided = m1 / 2.0;
    assert(m_divided.getValue() == 50.0);
    
    // Comparison
    assert(m1 > m2);
    assert(m2 < m1);
    assert(m1 != m2);
    
    std::cout << "  ✓ Passed" << std::endl;
}

// ===========================================
// Test 2: Compile-time dimension checking
// ===========================================
void test_dimension_checking() {
    std::cout << "Test 2: Compile-time dimension checking..." << std::endl;
    
    // These should compile:
    Mass<CGSUnit> m(100.0);
    Length<CGSUnit> l(50.0);
    
    // This should compile (same dimension):
    Mass<CGSUnit> m_sum = m + Mass<CGSUnit>(10.0);
    
    // This should NOT compile (uncomment to verify):
    // Mass<CGSUnit> error = m + l;  // Compile error: incompatible types!
    
    std::cout << "  ✓ Passed (compile-time check)" << std::endl;
}

// ===========================================
// Test 3: Unit conversions - CGS to SI
// ===========================================
void test_cgs_to_si_conversion() {
    std::cout << "Test 3: CGS to SI conversions..." << std::endl;
    
    // Mass: 1000 g = 1 kg
    Mass<CGSUnit> mass_cgs(1000.0);
    Mass<SIUnit> mass_si = convertUnit<SIUnit>(mass_cgs);
    assert(isClose(mass_si.getValue(), 1.0));
    
    // Length: 100 cm = 1 m
    Length<CGSUnit> length_cgs(100.0);
    Length<SIUnit> length_si = convertUnit<SIUnit>(length_cgs);
    assert(isClose(length_si.getValue(), 1.0));
    
    // Time: seconds are the same
    Time<CGSUnit> time_cgs(10.0);
    Time<SIUnit> time_si = convertUnit<SIUnit>(time_cgs);
    assert(isClose(time_si.getValue(), 10.0));
    
    // Velocity: 100 cm/s = 1 m/s
    Velocity<CGSUnit> vel_cgs(100.0);
    Velocity<SIUnit> vel_si = convertUnit<SIUnit>(vel_cgs);
    assert(isClose(vel_si.getValue(), 1.0));
    
    std::cout << "  ✓ Passed" << std::endl;
}

// ===========================================
// Test 4: Unit conversions - CGS to Galactic
// ===========================================
void test_cgs_to_galactic_conversion() {
    std::cout << "Test 4: CGS to Galactic conversions..." << std::endl;
    
    // Solar mass in CGS
    Mass<CGSUnit> solar_mass_cgs(CGSConstants::kSolarMass);
    Mass<GalacticUnit> solar_mass_gal = convertUnit<GalacticUnit>(solar_mass_cgs);
    assert(isClose(solar_mass_gal.getValue(), 1.0, 1e-6));
    
    // Parsec in CGS
    Length<CGSUnit> parsec_cgs(CGSConstants::kParsec);
    Length<GalacticUnit> parsec_gal = convertUnit<GalacticUnit>(parsec_cgs);
    assert(isClose(parsec_gal.getValue(), 1.0, 1e-6));
    
    // km/s in CGS
    Velocity<CGSUnit> kmps_cgs(CGSConstants::kKilometer);
    Velocity<GalacticUnit> kmps_gal = convertUnit<GalacticUnit>(kmps_cgs);
    assert(isClose(kmps_gal.getValue(), 1.0, 1e-6));
    
    std::cout << "  ✓ Passed" << std::endl;
}

// ===========================================
// Test 5: Round-trip conversions
// ===========================================
void test_roundtrip_conversions() {
    std::cout << "Test 5: Round-trip conversions..." << std::endl;
    
    // CGS -> SI -> CGS
    Mass<CGSUnit> original(123.456);
    Mass<SIUnit> si_version = convertUnit<SIUnit>(original);
    Mass<CGSUnit> back_to_cgs = convertUnit<CGSUnit>(si_version);
    assert(isClose(original.getValue(), back_to_cgs.getValue()));
    
    // CGS -> Galactic -> CGS
    Length<CGSUnit> length_original(CGSConstants::kKiloparsec * 5.0);
    Length<GalacticUnit> gal_version = convertUnit<GalacticUnit>(length_original);
    Length<CGSUnit> back_to_cgs_length = convertUnit<CGSUnit>(gal_version);
    assert(isClose(length_original.getValue(), back_to_cgs_length.getValue(), 1e-6));
    
    std::cout << "  ✓ Passed" << std::endl;
}

// ===========================================
// Test 6: Unit strings
// ===========================================
void test_unit_strings() {
    std::cout << "Test 6: Unit string representations..." << std::endl;
    
    // CGS
    assert(std::string(getUnitString<MassDim, CGSUnit>()) == "g");
    assert(std::string(getUnitString<LengthDim, CGSUnit>()) == "cm");
    assert(std::string(getUnitString<VelocityDim, CGSUnit>()) == "cm/s");
    
    // SI
    assert(std::string(getUnitString<MassDim, SIUnit>()) == "kg");
    assert(std::string(getUnitString<LengthDim, SIUnit>()) == "m");
    assert(std::string(getUnitString<VelocityDim, SIUnit>()) == "m/s");
    
    // Galactic
    assert(std::string(getUnitString<MassDim, GalacticUnit>()) == "M_sun");
    assert(std::string(getUnitString<LengthDim, GalacticUnit>()) == "pc");
    assert(std::string(getUnitString<VelocityDim, GalacticUnit>()) == "km/s");
    
    // Simulation
    assert(std::string(getUnitString<MassDim, SimulationUnit>()) == "sim");
    
    std::cout << "  ✓ Passed" << std::endl;
}

// ===========================================
// Test 7: CSV output creation
// ===========================================
void test_csv_output() {
    std::cout << "Test 7: CSV output..." << std::endl;
    
    // Create test particles
    std::vector<TestParticle> particles(10);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].id = i;
        particles[i].mass = 1.0 + i * 0.1;
        particles[i].pos = F64vec(i * 0.5, i * 0.3, i * 0.2);
        particles[i].vel = F64vec(i * 0.01, i * 0.02, i * 0.03);
        particles[i].dens = 1.0;
        particles[i].eng = 2.5;
        particles[i].pres = 1.5;
    }
    
    // Setup metadata
    OutputMetadata metadata;
    metadata.time = 0.5;
    metadata.total_particles = particles.size();
    metadata.mpi_rank = 0;
    metadata.mpi_size = 1;
    metadata.description = "Unit test";
    
    // Create CSV writer
    auto writer = createOutputWriter(
        OutputFormat::CSV,
        "test_output",
        metadata,
        CommInfo()
    );
    
    assert(writer != nullptr);
    
    // Write particles
    bool opened = writer->open();
    assert(opened);
    
    bool written = writer->writeParticles<TestParticle, SimulationUnit>(
        particles.data(), particles.size());
    assert(written);
    
    writer->close();
    
    std::cout << "  ✓ Passed (check test_output.csv)" << std::endl;
}

// ===========================================
// Test 8: Multiple unit systems in output
// ===========================================
void test_multiple_unit_outputs() {
    std::cout << "Test 8: Multiple unit system outputs..." << std::endl;
    
    // Create test particles
    std::vector<TestParticle> particles(5);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].id = i;
        particles[i].mass = CGSConstants::kSolarMass * (1.0 + i);
        particles[i].pos = F64vec(
            CGSConstants::kParsec * i,
            CGSConstants::kParsec * i * 0.5,
            CGSConstants::kParsec * i * 0.25
        );
        particles[i].vel = F64vec(
            CGSConstants::kKilometer * 100.0,
            CGSConstants::kKilometer * 50.0,
            CGSConstants::kKilometer * 25.0
        );
        particles[i].dens = 1.0e-24;  // g/cm^3
        particles[i].eng = 1.0e10;    // erg/g
        particles[i].pres = 1.0e-12;  // Ba
    }
    
    OutputMetadata metadata;
    metadata.time = 1.0;
    metadata.total_particles = particles.size();
    metadata.description = "Multi-unit test";
    
    // Output in CGS units
    {
        auto writer = createOutputWriter(
            OutputFormat::CSV, "test_cgs_output", metadata);
        writer->open();
        writer->writeParticles<TestParticle, CGSUnit>(
            particles.data(), particles.size());
        writer->close();
    }
    
    // Output in Galactic units
    {
        auto writer = createOutputWriter(
            OutputFormat::CSV, "test_galactic_output", metadata);
        writer->open();
        writer->writeParticles<TestParticle, GalacticUnit>(
            particles.data(), particles.size());
        writer->close();
    }
    
    // Output in SI units
    {
        auto writer = createOutputWriter(
            OutputFormat::CSV, "test_si_output", metadata);
        writer->open();
        writer->writeParticles<TestParticle, SIUnit>(
            particles.data(), particles.size());
        writer->close();
    }
    
    std::cout << "  ✓ Passed (check test_*_output.csv files)" << std::endl;
}

#ifdef PARTICLE_SIMULATOR_USE_HDF5
// ===========================================
// Test 9: HDF5 output (if available)
// ===========================================
void test_hdf5_output() {
    std::cout << "Test 9: HDF5 output..." << std::endl;
    
    // Create test particles
    std::vector<TestParticle> particles(100);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].id = i;
        particles[i].mass = 1.0 + i * 0.01;
        particles[i].pos = F64vec(
            std::sin(i * 0.1),
            std::cos(i * 0.1),
            i * 0.05
        );
        particles[i].vel = F64vec(
            std::cos(i * 0.1) * 0.1,
            -std::sin(i * 0.1) * 0.1,
            0.01
        );
        particles[i].dens = 1.0;
        particles[i].eng = 2.5;
        particles[i].pres = 1.5;
    }
    
    OutputMetadata metadata;
    metadata.time = 0.75;
    metadata.total_particles = particles.size();
    metadata.description = "HDF5 test";
    
    auto writer = createOutputWriter(
        OutputFormat::HDF5,
        "test_hdf5_output",
        metadata
    );
    
    assert(writer != nullptr);
    assert(writer->open());
    
    bool written = writer->writeParticles<TestParticle, GalacticUnit>(
        particles.data(), particles.size());
    assert(written);
    
    writer->close();
    
    std::cout << "  ✓ Passed (check test_hdf5_output.h5)" << std::endl;
}
#endif

// ===========================================
// Test 10: Physical constant verification
// ===========================================
void test_physical_constants() {
    std::cout << "Test 10: Physical constants..." << std::endl;
    
    // Check gravitational constant
    constexpr F64 G_cgs = CGSConstants::kGravitationalConstant;
    assert(isClose(G_cgs, 6.67430e-8, 1e-5));
    
    // Check solar mass
    constexpr F64 M_sun = CGSConstants::kSolarMass;
    assert(isClose(M_sun, 1.98892e33, 1e-5));
    
    // Check parsec
    constexpr F64 pc = CGSConstants::kParsec;
    assert(isClose(pc, 3.08567758149137e18, 1e-5));
    
    std::cout << "  ✓ Passed" << std::endl;
}

// ===========================================
// Main test runner
// ===========================================
int main(int argc, char* argv[]) {
    std::cout << "=========================================" << std::endl;
    std::cout << "FDPS Unit System & Output Tests" << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << std::endl;
    
    try {
        test_quantity_basic();
        test_dimension_checking();
        test_cgs_to_si_conversion();
        test_cgs_to_galactic_conversion();
        test_roundtrip_conversions();
        test_unit_strings();
        test_csv_output();
        test_multiple_unit_outputs();
        
#ifdef PARTICLE_SIMULATOR_USE_HDF5
        test_hdf5_output();
#else
        std::cout << "Test 9: HDF5 output... SKIPPED (HDF5 not available)" << std::endl;
#endif
        
        test_physical_constants();
        
        std::cout << std::endl;
        std::cout << "=========================================" << std::endl;
        std::cout << "All tests passed successfully! ✓" << std::endl;
        std::cout << "=========================================" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << std::endl;
        std::cerr << "TEST FAILED: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << std::endl;
        std::cerr << "TEST FAILED: Unknown exception" << std::endl;
        return 1;
    }
}
