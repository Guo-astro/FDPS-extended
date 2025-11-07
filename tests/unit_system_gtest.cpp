/*
 * Google Test suite for FDPS Unit System and Output Formats
 * 
 * This file tests:
 * 1. Compile-time dimension checking
 * 2. Unit conversions between different systems
 * 3. Output format functionality
 * 4. Extensibility of the system
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <vector>
#include <fstream>
#include "ps_defs.hpp"
#include "fdps_unit_system.hpp"
#include "fdps_output.hpp"

using namespace PS;
using namespace PS::UnitSystem;
using namespace PS::Output;
using ::testing::DoubleNear;
using ::testing::Not;

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

// Test fixture for common setup
class UnitSystemTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Setup code if needed
    }

    void TearDown() override {
        // Cleanup code if needed
    }

    // Helper function for floating point comparison
    bool isClose(F64 a, F64 b, F64 rel_tol = 1e-9, F64 abs_tol = 1e-12) {
        return std::abs(a - b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol);
    }
};

// ===========================================
// Test Suite 1: Basic Quantity Operations
// ===========================================

TEST_F(UnitSystemTest, QuantityAddition) {
    Mass<CGSUnit> m1(100.0);
    Mass<CGSUnit> m2(50.0);
    Mass<CGSUnit> m_sum = m1 + m2;
    
    EXPECT_DOUBLE_EQ(m_sum.getValue(), 150.0);
}

TEST_F(UnitSystemTest, QuantitySubtraction) {
    Mass<CGSUnit> m1(100.0);
    Mass<CGSUnit> m2(50.0);
    Mass<CGSUnit> m_diff = m1 - m2;
    
    EXPECT_DOUBLE_EQ(m_diff.getValue(), 50.0);
}

TEST_F(UnitSystemTest, ScalarMultiplication) {
    Mass<CGSUnit> m1(100.0);
    Mass<CGSUnit> m_scaled = m1 * 2.0;
    
    EXPECT_DOUBLE_EQ(m_scaled.getValue(), 200.0);
}

TEST_F(UnitSystemTest, ScalarDivision) {
    Mass<CGSUnit> m1(100.0);
    Mass<CGSUnit> m_divided = m1 / 2.0;
    
    EXPECT_DOUBLE_EQ(m_divided.getValue(), 50.0);
}

TEST_F(UnitSystemTest, ComparisonOperators) {
    Mass<CGSUnit> m1(100.0);
    Mass<CGSUnit> m2(50.0);
    
    EXPECT_GT(m1, m2);
    EXPECT_LT(m2, m1);
    EXPECT_NE(m1, m2);
    EXPECT_EQ(m1, Mass<CGSUnit>(100.0));
}

TEST_F(UnitSystemTest, CompoundAssignmentOperators) {
    Mass<CGSUnit> m(100.0);
    
    m += Mass<CGSUnit>(50.0);
    EXPECT_DOUBLE_EQ(m.getValue(), 150.0);
    
    m -= Mass<CGSUnit>(30.0);
    EXPECT_DOUBLE_EQ(m.getValue(), 120.0);
    
    m *= 2.0;
    EXPECT_DOUBLE_EQ(m.getValue(), 240.0);
    
    m /= 4.0;
    EXPECT_DOUBLE_EQ(m.getValue(), 60.0);
}

// ===========================================
// Test Suite 2: Unit Conversions - CGS to SI
// ===========================================

TEST_F(UnitSystemTest, MassConversionCGStoSI) {
    // 1000 g = 1 kg
    Mass<CGSUnit> mass_cgs(1000.0);
    Mass<SIUnit> mass_si = convertUnit<SIUnit>(mass_cgs);
    
    EXPECT_TRUE(isClose(mass_si.getValue(), 1.0));
}

TEST_F(UnitSystemTest, LengthConversionCGStoSI) {
    // 100 cm = 1 m
    Length<CGSUnit> length_cgs(100.0);
    Length<SIUnit> length_si = convertUnit<SIUnit>(length_cgs);
    
    EXPECT_TRUE(isClose(length_si.getValue(), 1.0));
}

TEST_F(UnitSystemTest, TimeConversionCGStoSI) {
    // seconds are the same
    Time<CGSUnit> time_cgs(10.0);
    Time<SIUnit> time_si = convertUnit<SIUnit>(time_cgs);
    
    EXPECT_DOUBLE_EQ(time_si.getValue(), 10.0);
}

TEST_F(UnitSystemTest, VelocityConversionCGStoSI) {
    // 100 cm/s = 1 m/s
    Velocity<CGSUnit> vel_cgs(100.0);
    Velocity<SIUnit> vel_si = convertUnit<SIUnit>(vel_cgs);
    
    EXPECT_TRUE(isClose(vel_si.getValue(), 1.0));
}

TEST_F(UnitSystemTest, EnergyConversionCGStoSI) {
    // 1 erg = 1e-7 J
    Energy<CGSUnit> energy_cgs(1.0);
    Energy<SIUnit> energy_si = convertUnit<SIUnit>(energy_cgs);
    
    EXPECT_TRUE(isClose(energy_si.getValue(), 1.0e-7, 1e-9));
}

TEST_F(UnitSystemTest, DensityConversionCGStoSI) {
    // 1 g/cm^3 = 1000 kg/m^3
    Density<CGSUnit> density_cgs(1.0);
    Density<SIUnit> density_si = convertUnit<SIUnit>(density_cgs);
    
    EXPECT_TRUE(isClose(density_si.getValue(), 1000.0));
}

// ===========================================
// Test Suite 3: Unit Conversions - CGS to Galactic
// ===========================================

TEST_F(UnitSystemTest, MassConversionCGStoGalactic) {
    // Solar mass in CGS should be 1.0 in Galactic units
    Mass<CGSUnit> solar_mass_cgs(CGSConstants::kSolarMass);
    Mass<GalacticUnit> solar_mass_gal = convertUnit<GalacticUnit>(solar_mass_cgs);
    
    EXPECT_TRUE(isClose(solar_mass_gal.getValue(), 1.0, 1e-6));
}

TEST_F(UnitSystemTest, LengthConversionCGStoGalactic) {
    // Parsec in CGS should be 1.0 in Galactic units
    Length<CGSUnit> parsec_cgs(CGSConstants::kParsec);
    Length<GalacticUnit> parsec_gal = convertUnit<GalacticUnit>(parsec_cgs);
    
    EXPECT_TRUE(isClose(parsec_gal.getValue(), 1.0, 1e-6));
}

TEST_F(UnitSystemTest, VelocityConversionCGStoGalactic) {
    // km/s in CGS should be 1.0 in Galactic units
    Velocity<CGSUnit> kmps_cgs(CGSConstants::kKilometer);
    Velocity<GalacticUnit> kmps_gal = convertUnit<GalacticUnit>(kmps_cgs);
    
    EXPECT_TRUE(isClose(kmps_gal.getValue(), 1.0, 1e-6));
}

// ===========================================
// Test Suite 4: Round-trip Conversions
// ===========================================

TEST_F(UnitSystemTest, RoundTripCGStoSItoCGS) {
    Mass<CGSUnit> original(123.456);
    Mass<SIUnit> si_version = convertUnit<SIUnit>(original);
    Mass<CGSUnit> back_to_cgs = convertUnit<CGSUnit>(si_version);
    
    EXPECT_TRUE(isClose(original.getValue(), back_to_cgs.getValue()));
}

TEST_F(UnitSystemTest, RoundTripCGStoGalactictoCGS) {
    Length<CGSUnit> length_original(CGSConstants::kKiloparsec * 5.0);
    Length<GalacticUnit> gal_version = convertUnit<GalacticUnit>(length_original);
    Length<CGSUnit> back_to_cgs_length = convertUnit<CGSUnit>(gal_version);
    
    EXPECT_TRUE(isClose(length_original.getValue(), back_to_cgs_length.getValue(), 1e-6));
}

TEST_F(UnitSystemTest, RoundTripSItoGalactictoSI) {
    Mass<SIUnit> mass_si(1000.0);  // 1000 kg
    Mass<CGSUnit> mass_cgs = convertUnit<CGSUnit>(mass_si);
    Mass<GalacticUnit> mass_gal = convertUnit<GalacticUnit>(mass_cgs);
    Mass<CGSUnit> back_cgs = convertUnit<CGSUnit>(mass_gal);
    Mass<SIUnit> back_si = convertUnit<SIUnit>(back_cgs);
    
    EXPECT_TRUE(isClose(mass_si.getValue(), back_si.getValue(), 1e-6));
}

// ===========================================
// Test Suite 5: Unit Strings
// ===========================================

TEST_F(UnitSystemTest, CGSUnitStrings) {
    const char* mass_str = getUnitString<MassDim, CGSUnit>();
    const char* length_str = getUnitString<LengthDim, CGSUnit>();
    const char* velocity_str = getUnitString<VelocityDim, CGSUnit>();
    const char* energy_str = getUnitString<EnergyDim, CGSUnit>();
    const char* density_str = getUnitString<DensityDim, CGSUnit>();
    
    EXPECT_STREQ(mass_str, "g");
    EXPECT_STREQ(length_str, "cm");
    EXPECT_STREQ(velocity_str, "cm/s");
    EXPECT_STREQ(energy_str, "erg");
    EXPECT_STREQ(density_str, "g/cm^3");
}

TEST_F(UnitSystemTest, SIUnitStrings) {
    const char* mass_str = getUnitString<MassDim, SIUnit>();
    const char* length_str = getUnitString<LengthDim, SIUnit>();
    const char* velocity_str = getUnitString<VelocityDim, SIUnit>();
    const char* energy_str = getUnitString<EnergyDim, SIUnit>();
    const char* density_str = getUnitString<DensityDim, SIUnit>();
    
    EXPECT_STREQ(mass_str, "kg");
    EXPECT_STREQ(length_str, "m");
    EXPECT_STREQ(velocity_str, "m/s");
    EXPECT_STREQ(energy_str, "J");
    EXPECT_STREQ(density_str, "kg/m^3");
}

TEST_F(UnitSystemTest, GalacticUnitStrings) {
    const char* mass_str = getUnitString<MassDim, GalacticUnit>();
    const char* length_str = getUnitString<LengthDim, GalacticUnit>();
    const char* velocity_str = getUnitString<VelocityDim, GalacticUnit>();
    
    EXPECT_STREQ(mass_str, "M_sun");
    EXPECT_STREQ(length_str, "pc");
    EXPECT_STREQ(velocity_str, "km/s");
}

TEST_F(UnitSystemTest, SimulationUnitStrings) {
    const char* mass_str = getUnitString<MassDim, SimulationUnit>();
    const char* length_str = getUnitString<LengthDim, SimulationUnit>();
    const char* velocity_str = getUnitString<VelocityDim, SimulationUnit>();
    
    EXPECT_STREQ(mass_str, "sim");
    EXPECT_STREQ(length_str, "sim");
    EXPECT_STREQ(velocity_str, "sim");
}

TEST_F(UnitSystemTest, UnitSystemNames) {
    const char* sim_name = getUnitSystemName<SimulationUnit>();
    const char* gal_name = getUnitSystemName<GalacticUnit>();
    const char* cgs_name = getUnitSystemName<CGSUnit>();
    const char* si_name = getUnitSystemName<SIUnit>();
    
    EXPECT_STREQ(sim_name, "Simulation");
    EXPECT_STREQ(gal_name, "Galactic");
    EXPECT_STREQ(cgs_name, "CGS");
    EXPECT_STREQ(si_name, "SI");
}

// ===========================================
// Test Suite 6: Physical Constants
// ===========================================

TEST_F(UnitSystemTest, GravitationalConstant) {
    constexpr F64 G_cgs = CGSConstants::kGravitationalConstant;
    EXPECT_TRUE(isClose(G_cgs, 6.67430e-8, 1e-5));
}

TEST_F(UnitSystemTest, SolarMass) {
    constexpr F64 M_sun = CGSConstants::kSolarMass;
    EXPECT_TRUE(isClose(M_sun, 1.98892e33, 1e-5));
}

TEST_F(UnitSystemTest, Parsec) {
    constexpr F64 pc = CGSConstants::kParsec;
    EXPECT_TRUE(isClose(pc, 3.08567758149137e18, 1e-5));
}

TEST_F(UnitSystemTest, Kilometer) {
    constexpr F64 km = CGSConstants::kKilometer;
    EXPECT_DOUBLE_EQ(km, 1.0e5);
}

// ===========================================
// Test Suite 7: CSV Output
// ===========================================

class OutputTest : public ::testing::Test {
protected:
    std::vector<TestParticle> createTestParticles(size_t count) {
        std::vector<TestParticle> particles(count);
        for (size_t i = 0; i < count; ++i) {
            particles[i].id = i;
            particles[i].mass = 1.0 + i * 0.1;
            particles[i].pos = F64vec(i * 0.5, i * 0.3, i * 0.2);
            particles[i].vel = F64vec(i * 0.01, i * 0.02, i * 0.03);
            particles[i].dens = 1.0;
            particles[i].eng = 2.5;
            particles[i].pres = 1.5;
        }
        return particles;
    }

    void TearDown() override {
        // Clean up test files
        std::remove("test_csv_basic.csv");
        std::remove("test_csv_sim.csv");
        std::remove("test_csv_cgs.csv");
        std::remove("test_csv_gal.csv");
    }
};

TEST_F(OutputTest, CSVWriterCreation) {
    OutputMetadata metadata;
    metadata.time = 0.5;
    metadata.total_particles = 10;
    metadata.description = "Test";
    
    auto writer = createOutputWriter(
        OutputFormat::CSV,
        "test_csv_basic",
        metadata,
        CommInfo()
    );
    
    ASSERT_NE(writer, nullptr);
}

TEST_F(OutputTest, CSVWriterOpenClose) {
    OutputMetadata metadata;
    metadata.time = 0.5;
    metadata.total_particles = 10;
    
    auto writer = createOutputWriter(
        OutputFormat::CSV,
        "test_csv_basic",
        metadata
    );
    
    ASSERT_TRUE(writer->open());
    EXPECT_TRUE(writer->isOpen());
    
    writer->close();
    EXPECT_FALSE(writer->isOpen());
}

TEST_F(OutputTest, CSVWriteParticles) {
    auto particles = createTestParticles(10);
    
    OutputMetadata metadata;
    metadata.time = 0.5;
    metadata.total_particles = particles.size();
    metadata.description = "Unit test";
    
    auto writer = createOutputWriter(
        OutputFormat::CSV,
        "test_csv_sim",
        metadata
    );
    
    ASSERT_TRUE(writer->open());
    bool write_success = writer->writeParticles<TestParticle, SimulationUnit>(
        particles.data(), particles.size());
    EXPECT_TRUE(write_success);
    writer->close();
    
    // Verify file exists
    std::ifstream file("test_csv_sim.csv");
    EXPECT_TRUE(file.good());
    file.close();
}

TEST_F(OutputTest, CSVMultipleUnitSystems) {
    auto particles = createTestParticles(5);
    
    OutputMetadata metadata;
    metadata.time = 1.0;
    metadata.total_particles = particles.size();
    metadata.description = "Multi-unit test";
    
    // Write in CGS
    {
        auto writer = createOutputWriter(
            OutputFormat::CSV, "test_csv_cgs", metadata);
        ASSERT_TRUE(writer->open());
        bool write_cgs = writer->writeParticles<TestParticle, CGSUnit>(
            particles.data(), particles.size());
        EXPECT_TRUE(write_cgs);
        writer->close();
    }
    
    // Write in Galactic
    {
        auto writer = createOutputWriter(
            OutputFormat::CSV, "test_csv_gal", metadata);
        ASSERT_TRUE(writer->open());
        bool write_gal = writer->writeParticles<TestParticle, GalacticUnit>(
            particles.data(), particles.size());
        EXPECT_TRUE(write_gal);
        writer->close();
    }
    
    // Verify both files exist
    std::ifstream file_cgs("test_csv_cgs.csv");
    EXPECT_TRUE(file_cgs.good());
    file_cgs.close();
    
    std::ifstream file_gal("test_csv_gal.csv");
    EXPECT_TRUE(file_gal.good());
    file_gal.close();
}

#ifdef PARTICLE_SIMULATOR_USE_HDF5
// ===========================================
// Test Suite 8: HDF5 Output (if available)
// ===========================================

class HDF5OutputTest : public ::testing::Test {
protected:
    std::vector<TestParticle> createTestParticles(size_t count) {
        std::vector<TestParticle> particles(count);
        for (size_t i = 0; i < count; ++i) {
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
        return particles;
    }

    void TearDown() override {
        std::remove("test_hdf5_output.h5");
    }
};

TEST_F(HDF5OutputTest, HDF5WriterCreation) {
    OutputMetadata metadata;
    metadata.time = 0.75;
    metadata.total_particles = 100;
    
    auto writer = createOutputWriter(
        OutputFormat::HDF5,
        "test_hdf5_output",
        metadata
    );
    
    ASSERT_NE(writer, nullptr);
}

TEST_F(HDF5OutputTest, HDF5WriteParticles) {
    auto particles = createTestParticles(100);
    
    OutputMetadata metadata;
    metadata.time = 0.75;
    metadata.total_particles = particles.size();
    metadata.description = "HDF5 test";
    
    auto writer = createOutputWriter(
        OutputFormat::HDF5,
        "test_hdf5_output",
        metadata
    );
    
    ASSERT_TRUE(writer->open());
    bool write_success = writer->writeParticles<TestParticle, GalacticUnit>(
        particles.data(), particles.size());
    EXPECT_TRUE(write_success);
    writer->close();
    
    // Verify file exists
    std::ifstream file("test_hdf5_output.h5", std::ios::binary);
    EXPECT_TRUE(file.good());
    file.close();
}
#endif

// Main function
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
