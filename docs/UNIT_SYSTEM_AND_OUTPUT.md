# FDPS Unit System and Output Format Extensions

## Overview

This implementation provides a **type-safe, compile-time checked unit system** and a **flexible, extensible output format system** for the Framework for Developing Particle Simulators (FDPS).

## Key Features

### 1. Type-Safe Unit System

- **Compile-time dimension checking**: Prevents mixing incompatible physical quantities
- **Multiple unit systems supported**:
  - **Simulation units** (dimensionless, G=M=R=1)
  - **Galactic units** (parsec, solar mass, km/s)
  - **CGS units** (cm, g, s)
  - **SI units** (m, kg, s)
- **Zero runtime overhead**: All conversions resolved at compile time when possible
- **Extensible**: Easy to add new unit systems

### 2. Flexible Output System

- **Multiple output formats**:
  - **CSV** with headers and unit metadata
  - **HDF5** with hierarchical structure and attributes
  - Extensible for future formats (VTK, NetCDF, etc.)
- **MPI-aware**: Automatically handles parallel output with rank suffixes
- **Strategy pattern**: Easy to switch between formats
- **Unit metadata**: All outputs include unit system information

## Architecture

### Unit System Design

```
Quantity<DimType, UnitSystem>
    ├── Compile-time dimension checking
    ├── Arithmetic operators (only for same dimensions)
    └── Unit conversion operators

UnitConverter<FromUnit, ToUnit>
    └── Template specializations for each conversion pair

Physical Dimensions:
    ├── MassDim
    ├── LengthDim
    ├── TimeDim
    ├── VelocityDim
    ├── EnergyDim
    ├── DensityDim
    ├── AccelerationDim
    └── PressureDim
```

### Output System Design

```
OutputWriter (Strategy Pattern)
    ├── CSVWriter
    ├── HDF5Writer (optional, requires HDF5 library)
    └── Future writers...

OutputMetadata
    ├── time
    ├── total_particles
    ├── unit_system_name
    └── description
```

## Usage Examples

### Basic Unit System Usage

```cpp
#include "fdps_unit_system.hpp"

using namespace PS::UnitSystem;

// Create quantities in CGS units
Mass<CGSUnit> mass_cgs(1.98892e33);  // Solar mass in grams
Length<CGSUnit> length_cgs(3.086e18); // 1 parsec in cm

// Convert to Galactic units
auto mass_gal = convertUnit<GalacticUnit>(mass_cgs);
auto length_gal = convertUnit<GalacticUnit>(length_cgs);

// These operations are type-safe:
Mass<CGSUnit> total_mass = mass_cgs + Mass<CGSUnit>(1.0e30); // OK
// Mass<CGSUnit> error = mass_cgs + length_cgs; // Compile error!

// Get unit string for output
const char* unit_str = getUnitString<MassDim, GalacticUnit>(); // "M_sun"
```

### Using CSV Output

```cpp
#include "fdps_output.hpp"

using namespace PS::Output;
using namespace PS::UnitSystem;

// Setup metadata
OutputMetadata metadata;
metadata.time = 0.5;
metadata.total_particles = 10000;
metadata.description = "SPH shock tube test";

// Create CSV writer
auto writer = createOutputWriter(
    OutputFormat::CSV,
    "result/snapshot_0001",
    metadata,
    PS::CommInfo()
);

// Open and write
if (writer->open()) {
    const FP* particles = system.getParticlePointer();
    S64 num_local = system.getNumberOfParticleLocal();
    
    // Write with specific unit system
    writer->writeParticles<FP, SimulationUnit>(particles, num_local);
    writer->close();
}
```

### Using HDF5 Output

```cpp
#ifdef PARTICLE_SIMULATOR_USE_HDF5
// Create HDF5 writer (same interface!)
auto hdf5_writer = createOutputWriter(
    OutputFormat::HDF5,
    "result/snapshot_0001",
    metadata,
    PS::CommInfo()
);

if (hdf5_writer->open()) {
    hdf5_writer->writeParticles<FP, GalacticUnit>(particles, num_local);
    hdf5_writer->close();
}
#endif
```

### Complete Example: Multiple Outputs

```cpp
// Output the same data in multiple formats and unit systems
template<typename UnitSystem>
void outputSnapshot(const PS::ParticleSystem<FP>& system,
                   PS::Output::OutputFormat format,
                   const char* base_name,
                   PS::F64 time)
{
    using namespace PS::Output;
    
    OutputMetadata meta;
    meta.time = time;
    meta.total_particles = system.getNumberOfParticleGlobal();
    meta.mpi_rank = PS::Comm::getRank();
    meta.mpi_size = PS::Comm::getNumberOfProc();
    
    auto writer = createOutputWriter(format, base_name, meta);
    if (writer && writer->open()) {
        writer->writeParticles<FP, UnitSystem>(
            system.getParticlePointer(),
            system.getNumberOfParticleLocal()
        );
        writer->close();
    }
}

// In main simulation loop:
if (step % OUTPUT_INTERVAL == 0) {
    char filename[256];
    
    // CSV output in simulation units
    sprintf(filename, "result/snap_%04d_sim", step);
    outputSnapshot<SimulationUnit>(sph_system, OutputFormat::CSV, filename, time);
    
    // CSV output in galactic units
    sprintf(filename, "result/snap_%04d_gal", step);
    outputSnapshot<GalacticUnit>(sph_system, OutputFormat::CSV, filename, time);
    
    // HDF5 output in CGS units
    sprintf(filename, "result/snap_%04d_cgs", step);
    outputSnapshot<CGSUnit>(sph_system, OutputFormat::HDF5, filename, time);
}
```

## CSV Output Format

CSV files include:
- Metadata as comment lines (prefixed with `#`)
- Column headers
- Unit information for each column
- High-precision scientific notation (15 digits)

Example CSV output:
```csv
# Particle Simulation Output (CSV Format)
# Time: 5.000000000000000e-01
# Total Particles: 16384
# Unit System: Galactic
# MPI Rank: 0 of 4
#
id,mass,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,dens,eng,pres
# Units: none,M_sun,pc,pc,pc,km/s,km/s,km/s,M_sun/pc^3,M_sun*(km/s)^2,M_sun/(pc*(km/s)^2)
0,7.500000000000000e-01,1.234567890123456e-02,...
```

## HDF5 Output Format

HDF5 files include:
- `/Particles/` group containing all particle datasets
- Separate datasets for each field (Mass, PositionX, VelocityX, etc.)
- Unit information stored as dataset attributes
- Metadata stored as file attributes
- Efficient binary storage suitable for large simulations

Structure:
```
file.h5
├── Attributes: Time, TotalParticles, UnitSystem, MPIRank, MPISize
└── /Particles/
    ├── ID (int64 dataset)
    ├── Mass (float64 dataset, Units="M_sun")
    ├── PositionX (float64 dataset, Units="pc")
    ├── PositionY (float64 dataset, Units="pc")
    ├── PositionZ (float64 dataset, Units="pc")
    ├── VelocityX (float64 dataset, Units="km/s")
    ├── VelocityY (float64 dataset, Units="km/s")
    ├── VelocityZ (float64 dataset, Units="km/s")
    ├── Density (float64 dataset, Units="M_sun/pc^3")
    ├── Energy (float64 dataset, Units="M_sun*(km/s)^2")
    └── Pressure (float64 dataset, Units="M_sun/(pc*(km/s)^2)")
```

## Adding New Unit Systems

To add a new unit system (e.g., atomic units):

1. **Define unit system tag**:
```cpp
struct AtomicUnit {};
```

2. **Specify conversion factors**:
```cpp
namespace ConversionFactors {
    constexpr F64 kAtomicUnit_Mass = 9.10938e-28;  // electron mass in g
    constexpr F64 kAtomicUnit_Length = 5.29177e-9;  // Bohr radius in cm
    // ... etc
}
```

3. **Create UnitConverter specializations**:
```cpp
template<>
class UnitConverter<CGSUnit, AtomicUnit> {
public:
    static constexpr Mass<AtomicUnit> convert(const Mass<CGSUnit>& mass) noexcept {
        return Mass<AtomicUnit>(mass.getValue() / ConversionFactors::kAtomicUnit_Mass);
    }
    // ... other conversions
};
```

4. **Add unit strings**:
```cpp
template<> inline const char* getUnitString<MassDim, AtomicUnit>() noexcept { 
    return "m_e"; 
}
template<> inline const char* getUnitString<LengthDim, AtomicUnit>() noexcept { 
    return "a_0"; 
}
// ... etc
```

## Adding New Output Formats

To add a new output format (e.g., VTK):

1. **Create writer class**:
```cpp
class VTKWriter : public OutputWriter {
private:
    // VTK-specific members
    
public:
    VTKWriter(const char* filename, const OutputMetadata& metadata, 
              const CommInfo& comm_info)
        : OutputWriter(filename, metadata, comm_info) {}
    
    bool open() override {
        // VTK file opening logic
    }
    
    void close() override {
        // VTK file closing logic
    }
    
    bool isOpen() const override {
        // Check if file is open
    }
    
protected:
    bool writeHeader() override {
        // Write VTK header
    }
    
    bool writeFooter() override {
        // Write VTK footer
    }
    
    template<typename ParticleType, typename UnitSystem>
    friend bool OutputWriter::writeParticleData(...);
    
private:
    template<typename ParticleType, typename UnitSystem>
    bool writeParticleDataImpl(const ParticleType* particles, S64 num_particles) {
        // VTK-specific particle data writing
    }
};
```

2. **Update factory function**:
```cpp
inline std::unique_ptr<OutputWriter> createOutputWriter(...) {
    switch (format) {
        case OutputFormat::CSV:
            return std::make_unique<CSVWriter>(...);
        case OutputFormat::VTK:
            return std::make_unique<VTKWriter>(...);
        // ...
    }
}
```

3. **Add to OutputFormat enum**:
```cpp
enum class OutputFormat {
    CSV,
    HDF5,
    VTK,  // New format
    // ...
};
```

## Compilation

### Basic compilation (CSV only):
```bash
cd sample/c++/sph
g++ -std=c++17 -O3 -I../../../src main_with_units.cpp -o sph_with_units
```

### With HDF5 support:
```bash
g++ -std=c++17 -O3 -I../../../src -DPARTICLE_SIMULATOR_USE_HDF5 \
    main_with_units.cpp -lhdf5 -o sph_with_units
```

### With MPI:
```bash
mpic++ -std=c++17 -O3 -I../../../src -DPARTICLE_SIMULATOR_MPI_PARALLEL \
       main_with_units.cpp -o sph_with_units
```

### With MPI and HDF5:
```bash
mpic++ -std=c++17 -O3 -I../../../src \
       -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_USE_HDF5 \
       main_with_units.cpp -lhdf5 -o sph_with_units
```

## Design Principles

### Type Safety
- **Compile-time checking**: Dimension mismatches caught at compile time
- **No implicit conversions**: Must explicitly convert between unit systems
- **Strong typing**: Each quantity has both dimension and unit system information

### Zero-Cost Abstraction
- **Template metaprogramming**: No runtime overhead for unit checking
- **Inline functions**: Conversions optimized away when possible
- **constexpr**: Compile-time evaluation where applicable

### Extensibility
- **Open-Closed Principle**: Easy to add new units/formats without modifying existing code
- **Strategy Pattern**: Output formats selected at runtime
- **Template-based**: Type-safe and efficient

### Modularity
- **Single Responsibility**: Each class has one clear purpose
- **Separation of Concerns**: Units, output, and simulation logic are independent
- **Header-only**: Easy to integrate (for most components)

## Performance Considerations

- **Unit system**: Zero runtime overhead when using template-based approach
- **CSV output**: Text-based, slower but human-readable and portable
- **HDF5 output**: Binary, fast, and suitable for large datasets
- **MPI I/O**: Each rank writes to separate file (simplicity over performance)
  - For production, consider parallel HDF5 with collective I/O

## Future Extensions

Potential additions:
1. **More unit systems**: Planck units, natural units, etc.
2. **More output formats**: VTK, NetCDF, XDMF, etc.
3. **Parallel HDF5**: Collective MPI I/O to single file
4. **Compression**: gzip for CSV, compression filters for HDF5
5. **Checkpoint/restart**: Full state serialization
6. **Visualization helpers**: Python scripts to read outputs
7. **Unit conversion utilities**: Command-line tools for post-processing

## References

- FDPS Documentation: https://github.com/FDPS/FDPS
- HDF5 Documentation: https://www.hdfgroup.org/
- Physical Constants: CODATA 2018 values
- Unit Systems: IAU standard astronomical constants

## Author

Implementation following FDPS coding standards and best practices for C++ template metaprogramming.
