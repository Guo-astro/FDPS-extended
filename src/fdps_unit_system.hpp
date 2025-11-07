#pragma once

#include <type_traits>
#include <cmath>
#include <ratio>
#include "ps_defs.hpp"

namespace ParticleSimulator {
namespace UnitSystem {

// Dimension tags for compile-time dimension checking
struct MassDim {};
struct LengthDim {};
struct TimeDim {};
struct VelocityDim {};
struct EnergyDim {};
struct DensityDim {};
struct AccelerationDim {};
struct PressureDim {};

// Unit system tags
struct SimulationUnit {};
struct GalacticUnit {};
struct CGSUnit {};
struct SIUnit {};

// Compile-time dimensional analysis
template<int M, int L, int T>
struct Dimension {
    static constexpr int mass = M;
    static constexpr int length = L;
    static constexpr int time = T;
};

// Physical quantity with compile-time dimension checking
template<typename DimType, typename UnitSystem>
class Quantity {
private:
    F64 value_;

public:
    using dimension_type = DimType;
    using unit_system = UnitSystem;

    constexpr explicit Quantity(F64 val = 0.0) noexcept : value_(val) {}

    constexpr F64 getValue() const noexcept { return value_; }
    constexpr void setValue(F64 val) noexcept { value_ = val; }

    // Arithmetic operators with same dimension and unit system
    constexpr Quantity operator+(const Quantity& other) const noexcept {
        return Quantity(value_ + other.value_);
    }

    constexpr Quantity operator-(const Quantity& other) const noexcept {
        return Quantity(value_ - other.value_);
    }

    constexpr Quantity& operator+=(const Quantity& other) noexcept {
        value_ += other.value_;
        return *this;
    }

    constexpr Quantity& operator-=(const Quantity& other) noexcept {
        value_ -= other.value_;
        return *this;
    }

    // Scalar multiplication
    constexpr Quantity operator*(F64 scalar) const noexcept {
        return Quantity(value_ * scalar);
    }

    constexpr Quantity operator/(F64 scalar) const noexcept {
        return Quantity(value_ / scalar);
    }

    constexpr Quantity& operator*=(F64 scalar) noexcept {
        value_ *= scalar;
        return *this;
    }

    constexpr Quantity& operator/=(F64 scalar) noexcept {
        value_ /= scalar;
        return *this;
    }

    // Comparison operators
    constexpr bool operator==(const Quantity& other) const noexcept {
        return value_ == other.value_;
    }

    constexpr bool operator!=(const Quantity& other) const noexcept {
        return value_ != other.value_;
    }

    constexpr bool operator<(const Quantity& other) const noexcept {
        return value_ < other.value_;
    }

    constexpr bool operator>(const Quantity& other) const noexcept {
        return value_ > other.value_;
    }

    constexpr bool operator<=(const Quantity& other) const noexcept {
        return value_ <= other.value_;
    }

    constexpr bool operator>=(const Quantity& other) const noexcept {
        return value_ >= other.value_;
    }

    // Unary operators
    constexpr Quantity operator-() const noexcept {
        return Quantity(-value_);
    }
};

// Scalar multiplication (scalar * quantity)
template<typename DimType, typename UnitSystem>
constexpr Quantity<DimType, UnitSystem> operator*(F64 scalar, const Quantity<DimType, UnitSystem>& q) noexcept {
    return q * scalar;
}

// Type aliases for common quantities in each unit system
template<typename UnitSystem>
using Mass = Quantity<MassDim, UnitSystem>;

template<typename UnitSystem>
using Length = Quantity<LengthDim, UnitSystem>;

template<typename UnitSystem>
using Time = Quantity<TimeDim, UnitSystem>;

template<typename UnitSystem>
using Velocity = Quantity<VelocityDim, UnitSystem>;

template<typename UnitSystem>
using Energy = Quantity<EnergyDim, UnitSystem>;

template<typename UnitSystem>
using Density = Quantity<DensityDim, UnitSystem>;

template<typename UnitSystem>
using Acceleration = Quantity<AccelerationDim, UnitSystem>;

template<typename UnitSystem>
using Pressure = Quantity<PressureDim, UnitSystem>;

// Physical constants in CGS units (base)
namespace CGSConstants {
    constexpr F64 kGravitationalConstant = 6.67430e-8;    // cm^3/(g·s^2)
    constexpr F64 kSolarMass = 1.98892e33;                // g
    constexpr F64 kParsec = 3.08567758149137e18;          // cm
    constexpr F64 kKilometer = 1.0e5;                     // cm
    constexpr F64 kSecond = 1.0;                          // s
    constexpr F64 kYear = 3.15576e7;                      // s
    constexpr F64 kMegayear = 1.0e6 * kYear;              // s
    constexpr F64 kGigayear = 1.0e9 * kYear;              // s
    constexpr F64 kKiloparsec = 1.0e3 * kParsec;          // cm
    constexpr F64 kMegaparsec = 1.0e6 * kParsec;          // cm
}

// Conversion factors between unit systems
namespace ConversionFactors {
    // From CGS to other systems
    constexpr F64 kCGSToSI_Mass = 1.0e-3;                 // g -> kg
    constexpr F64 kCGSToSI_Length = 1.0e-2;               // cm -> m
    constexpr F64 kCGSToSI_Time = 1.0;                    // s -> s

    // Galactic units: parsec, solar mass, km/s
    constexpr F64 kGalacticUnit_Mass = CGSConstants::kSolarMass;
    constexpr F64 kGalacticUnit_Length = CGSConstants::kParsec;
    constexpr F64 kGalacticUnit_Velocity = CGSConstants::kKilometer;
    constexpr F64 kGalacticUnit_Time = kGalacticUnit_Length / kGalacticUnit_Velocity;
}

// Unit converter class
template<typename FromUnit, typename ToUnit>
class UnitConverter {
public:
    template<typename DimType>
    static constexpr Quantity<DimType, ToUnit> convert(const Quantity<DimType, DimType>& from) noexcept;
};

// Specialization: CGS to SI
template<>
class UnitConverter<CGSUnit, SIUnit> {
public:
    static constexpr Mass<SIUnit> convert(const Mass<CGSUnit>& mass) noexcept {
        return Mass<SIUnit>(mass.getValue() * ConversionFactors::kCGSToSI_Mass);
    }

    static constexpr Length<SIUnit> convert(const Length<CGSUnit>& length) noexcept {
        return Length<SIUnit>(length.getValue() * ConversionFactors::kCGSToSI_Length);
    }

    static constexpr Time<SIUnit> convert(const Time<CGSUnit>& time) noexcept {
        return Time<SIUnit>(time.getValue() * ConversionFactors::kCGSToSI_Time);
    }

    static constexpr Velocity<SIUnit> convert(const Velocity<CGSUnit>& vel) noexcept {
        return Velocity<SIUnit>(vel.getValue() * ConversionFactors::kCGSToSI_Length / 
                                                  ConversionFactors::kCGSToSI_Time);
    }

    static constexpr Energy<SIUnit> convert(const Energy<CGSUnit>& energy) noexcept {
        return Energy<SIUnit>(energy.getValue() * ConversionFactors::kCGSToSI_Mass * 
                                                   ConversionFactors::kCGSToSI_Length * 
                                                   ConversionFactors::kCGSToSI_Length /
                                                   (ConversionFactors::kCGSToSI_Time * 
                                                    ConversionFactors::kCGSToSI_Time));
    }

    static constexpr Density<SIUnit> convert(const Density<CGSUnit>& density) noexcept {
        return Density<SIUnit>(density.getValue() * ConversionFactors::kCGSToSI_Mass /
                                                     (ConversionFactors::kCGSToSI_Length *
                                                      ConversionFactors::kCGSToSI_Length *
                                                      ConversionFactors::kCGSToSI_Length));
    }

    static constexpr Acceleration<SIUnit> convert(const Acceleration<CGSUnit>& acc) noexcept {
        return Acceleration<SIUnit>(acc.getValue() * ConversionFactors::kCGSToSI_Length /
                                                      (ConversionFactors::kCGSToSI_Time *
                                                       ConversionFactors::kCGSToSI_Time));
    }

    static constexpr Pressure<SIUnit> convert(const Pressure<CGSUnit>& pressure) noexcept {
        return Pressure<SIUnit>(pressure.getValue() * ConversionFactors::kCGSToSI_Mass /
                                                       (ConversionFactors::kCGSToSI_Length *
                                                        ConversionFactors::kCGSToSI_Time *
                                                        ConversionFactors::kCGSToSI_Time));
    }
};

// Specialization: CGS to Galactic
template<>
class UnitConverter<CGSUnit, GalacticUnit> {
public:
    static constexpr Mass<GalacticUnit> convert(const Mass<CGSUnit>& mass) noexcept {
        return Mass<GalacticUnit>(mass.getValue() / ConversionFactors::kGalacticUnit_Mass);
    }

    static constexpr Length<GalacticUnit> convert(const Length<CGSUnit>& length) noexcept {
        return Length<GalacticUnit>(length.getValue() / ConversionFactors::kGalacticUnit_Length);
    }

    static constexpr Time<GalacticUnit> convert(const Time<CGSUnit>& time) noexcept {
        return Time<GalacticUnit>(time.getValue() / ConversionFactors::kGalacticUnit_Time);
    }

    static constexpr Velocity<GalacticUnit> convert(const Velocity<CGSUnit>& vel) noexcept {
        return Velocity<GalacticUnit>(vel.getValue() / ConversionFactors::kGalacticUnit_Velocity);
    }

    static constexpr Energy<GalacticUnit> convert(const Energy<CGSUnit>& energy) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Mass * 
                           ConversionFactors::kGalacticUnit_Velocity * 
                           ConversionFactors::kGalacticUnit_Velocity;
        return Energy<GalacticUnit>(energy.getValue() / factor);
    }

    static constexpr Density<GalacticUnit> convert(const Density<CGSUnit>& density) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Mass / 
                           (ConversionFactors::kGalacticUnit_Length *
                            ConversionFactors::kGalacticUnit_Length *
                            ConversionFactors::kGalacticUnit_Length);
        return Density<GalacticUnit>(density.getValue() / factor);
    }

    static constexpr Acceleration<GalacticUnit> convert(const Acceleration<CGSUnit>& acc) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Velocity / 
                           ConversionFactors::kGalacticUnit_Time;
        return Acceleration<GalacticUnit>(acc.getValue() / factor);
    }

    static constexpr Pressure<GalacticUnit> convert(const Pressure<CGSUnit>& pressure) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Mass / 
                           (ConversionFactors::kGalacticUnit_Length * 
                            ConversionFactors::kGalacticUnit_Time * 
                            ConversionFactors::kGalacticUnit_Time);
        return Pressure<GalacticUnit>(pressure.getValue() / factor);
    }
};

// Specialization: SI to CGS
template<>
class UnitConverter<SIUnit, CGSUnit> {
public:
    static constexpr Mass<CGSUnit> convert(const Mass<SIUnit>& mass) noexcept {
        return Mass<CGSUnit>(mass.getValue() / ConversionFactors::kCGSToSI_Mass);
    }

    static constexpr Length<CGSUnit> convert(const Length<SIUnit>& length) noexcept {
        return Length<CGSUnit>(length.getValue() / ConversionFactors::kCGSToSI_Length);
    }

    static constexpr Time<CGSUnit> convert(const Time<SIUnit>& time) noexcept {
        return Time<CGSUnit>(time.getValue() / ConversionFactors::kCGSToSI_Time);
    }

    static constexpr Velocity<CGSUnit> convert(const Velocity<SIUnit>& vel) noexcept {
        return Velocity<CGSUnit>(vel.getValue() / (ConversionFactors::kCGSToSI_Length / 
                                                    ConversionFactors::kCGSToSI_Time));
    }

    static constexpr Energy<CGSUnit> convert(const Energy<SIUnit>& energy) noexcept {
        const F64 factor = ConversionFactors::kCGSToSI_Mass * 
                           ConversionFactors::kCGSToSI_Length * 
                           ConversionFactors::kCGSToSI_Length /
                           (ConversionFactors::kCGSToSI_Time * 
                            ConversionFactors::kCGSToSI_Time);
        return Energy<CGSUnit>(energy.getValue() / factor);
    }

    static constexpr Density<CGSUnit> convert(const Density<SIUnit>& density) noexcept {
        const F64 factor = ConversionFactors::kCGSToSI_Mass /
                           (ConversionFactors::kCGSToSI_Length *
                            ConversionFactors::kCGSToSI_Length *
                            ConversionFactors::kCGSToSI_Length);
        return Density<CGSUnit>(density.getValue() / factor);
    }

    static constexpr Acceleration<CGSUnit> convert(const Acceleration<SIUnit>& acc) noexcept {
        const F64 factor = ConversionFactors::kCGSToSI_Length /
                           (ConversionFactors::kCGSToSI_Time *
                            ConversionFactors::kCGSToSI_Time);
        return Acceleration<CGSUnit>(acc.getValue() / factor);
    }

    static constexpr Pressure<CGSUnit> convert(const Pressure<SIUnit>& pressure) noexcept {
        const F64 factor = ConversionFactors::kCGSToSI_Mass /
                           (ConversionFactors::kCGSToSI_Length *
                            ConversionFactors::kCGSToSI_Time *
                            ConversionFactors::kCGSToSI_Time);
        return Pressure<CGSUnit>(pressure.getValue() / factor);
    }
};

// Specialization: Galactic to CGS
template<>
class UnitConverter<GalacticUnit, CGSUnit> {
public:
    static constexpr Mass<CGSUnit> convert(const Mass<GalacticUnit>& mass) noexcept {
        return Mass<CGSUnit>(mass.getValue() * ConversionFactors::kGalacticUnit_Mass);
    }

    static constexpr Length<CGSUnit> convert(const Length<GalacticUnit>& length) noexcept {
        return Length<CGSUnit>(length.getValue() * ConversionFactors::kGalacticUnit_Length);
    }

    static constexpr Time<CGSUnit> convert(const Time<GalacticUnit>& time) noexcept {
        return Time<CGSUnit>(time.getValue() * ConversionFactors::kGalacticUnit_Time);
    }

    static constexpr Velocity<CGSUnit> convert(const Velocity<GalacticUnit>& vel) noexcept {
        return Velocity<CGSUnit>(vel.getValue() * ConversionFactors::kGalacticUnit_Velocity);
    }

    static constexpr Energy<CGSUnit> convert(const Energy<GalacticUnit>& energy) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Mass * 
                           ConversionFactors::kGalacticUnit_Velocity * 
                           ConversionFactors::kGalacticUnit_Velocity;
        return Energy<CGSUnit>(energy.getValue() * factor);
    }

    static constexpr Density<CGSUnit> convert(const Density<GalacticUnit>& density) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Mass / 
                           (ConversionFactors::kGalacticUnit_Length *
                            ConversionFactors::kGalacticUnit_Length *
                            ConversionFactors::kGalacticUnit_Length);
        return Density<CGSUnit>(density.getValue() * factor);
    }

    static constexpr Acceleration<CGSUnit> convert(const Acceleration<GalacticUnit>& acc) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Velocity / 
                           ConversionFactors::kGalacticUnit_Time;
        return Acceleration<CGSUnit>(acc.getValue() * factor);
    }

    static constexpr Pressure<CGSUnit> convert(const Pressure<GalacticUnit>& pressure) noexcept {
        const F64 factor = ConversionFactors::kGalacticUnit_Mass / 
                           (ConversionFactors::kGalacticUnit_Length * 
                            ConversionFactors::kGalacticUnit_Time * 
                            ConversionFactors::kGalacticUnit_Time);
        return Pressure<CGSUnit>(pressure.getValue() * factor);
    }
};

// Helper function for convenient conversion
template<typename ToUnit, typename DimType, typename FromUnit>
constexpr auto convertUnit(const Quantity<DimType, FromUnit>& quantity) noexcept 
    -> Quantity<DimType, ToUnit> {
    return UnitConverter<FromUnit, ToUnit>::convert(quantity);
}

// Unit system name strings for output
template<typename UnitSystem>
inline const char* getUnitSystemName() noexcept;

template<>
inline const char* getUnitSystemName<SimulationUnit>() noexcept {
    return "Simulation";
}

template<>
inline const char* getUnitSystemName<GalacticUnit>() noexcept {
    return "Galactic";
}

template<>
inline const char* getUnitSystemName<CGSUnit>() noexcept {
    return "CGS";
}

template<>
inline const char* getUnitSystemName<SIUnit>() noexcept {
    return "SI";
}

// Unit strings for each dimension and unit system
template<typename DimType, typename UnitSystem>
inline const char* getUnitString() noexcept;

// CGS units
template<> inline const char* getUnitString<MassDim, CGSUnit>() noexcept { return "g"; }
template<> inline const char* getUnitString<LengthDim, CGSUnit>() noexcept { return "cm"; }
template<> inline const char* getUnitString<TimeDim, CGSUnit>() noexcept { return "s"; }
template<> inline const char* getUnitString<VelocityDim, CGSUnit>() noexcept { return "cm/s"; }
template<> inline const char* getUnitString<EnergyDim, CGSUnit>() noexcept { return "erg"; }
template<> inline const char* getUnitString<DensityDim, CGSUnit>() noexcept { return "g/cm^3"; }
template<> inline const char* getUnitString<AccelerationDim, CGSUnit>() noexcept { return "cm/s^2"; }
template<> inline const char* getUnitString<PressureDim, CGSUnit>() noexcept { return "Ba"; }

// SI units
template<> inline const char* getUnitString<MassDim, SIUnit>() noexcept { return "kg"; }
template<> inline const char* getUnitString<LengthDim, SIUnit>() noexcept { return "m"; }
template<> inline const char* getUnitString<TimeDim, SIUnit>() noexcept { return "s"; }
template<> inline const char* getUnitString<VelocityDim, SIUnit>() noexcept { return "m/s"; }
template<> inline const char* getUnitString<EnergyDim, SIUnit>() noexcept { return "J"; }
template<> inline const char* getUnitString<DensityDim, SIUnit>() noexcept { return "kg/m^3"; }
template<> inline const char* getUnitString<AccelerationDim, SIUnit>() noexcept { return "m/s^2"; }
template<> inline const char* getUnitString<PressureDim, SIUnit>() noexcept { return "Pa"; }

// Galactic units
template<> inline const char* getUnitString<MassDim, GalacticUnit>() noexcept { return "M_sun"; }
template<> inline const char* getUnitString<LengthDim, GalacticUnit>() noexcept { return "pc"; }
template<> inline const char* getUnitString<TimeDim, GalacticUnit>() noexcept { return "pc/(km/s)"; }
template<> inline const char* getUnitString<VelocityDim, GalacticUnit>() noexcept { return "km/s"; }
template<> inline const char* getUnitString<EnergyDim, GalacticUnit>() noexcept { return "M_sun*(km/s)^2"; }
template<> inline const char* getUnitString<DensityDim, GalacticUnit>() noexcept { return "M_sun/pc^3"; }
template<> inline const char* getUnitString<AccelerationDim, GalacticUnit>() noexcept { return "(km/s)^2/pc"; }
template<> inline const char* getUnitString<PressureDim, GalacticUnit>() noexcept { return "M_sun/(pc*(km/s)^2)"; }

// Simulation units (dimensionless)
template<> inline const char* getUnitString<MassDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<LengthDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<TimeDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<VelocityDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<EnergyDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<DensityDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<AccelerationDim, SimulationUnit>() noexcept { return "sim"; }
template<> inline const char* getUnitString<PressureDim, SimulationUnit>() noexcept { return "sim"; }

} // namespace UnitSystem
} // namespace ParticleSimulator
