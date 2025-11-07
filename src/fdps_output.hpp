#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "ps_defs.hpp"
#include "fdps_unit_system.hpp"

namespace ParticleSimulator {
namespace Output {

// Output format enum for easy selection
enum class OutputFormat {
    CSV,
    HDF5,
    Binary,
    ASCII
};

// Metadata structure for output files
struct OutputMetadata {
    F64 time;
    S64 total_particles;
    S32 mpi_rank;
    S32 mpi_size;
    std::string unit_system_name;
    std::string time_unit;
    std::string description;
};

// Abstract base class for output writers (Strategy pattern)
class OutputWriter {
protected:
    std::string filename_;
    OutputMetadata metadata_;
    CommInfo comm_info_;

public:
    OutputWriter(const char* filename, const OutputMetadata& metadata, const CommInfo& comm_info)
        : filename_(filename), metadata_(metadata), comm_info_(comm_info) {}

    virtual ~OutputWriter() = default;

    // Delete copy operations to prevent issues with file handles
    OutputWriter(const OutputWriter&) = delete;
    OutputWriter& operator=(const OutputWriter&) = delete;

    // Move operations
    OutputWriter(OutputWriter&&) = default;
    OutputWriter& operator=(OutputWriter&&) = default;

    // Pure virtual methods for derived classes to implement
    virtual bool open() = 0;
    virtual void close() = 0;
    virtual bool isOpen() const = 0;

    // Write particle data (template method pattern)
    template<typename ParticleType, typename UnitSystem>
    bool writeParticles(const ParticleType* particles, S64 num_particles);

protected:
    // Hook methods for derived classes
    virtual bool writeHeader() = 0;
    virtual bool writeFooter() = 0;

    template<typename ParticleType, typename UnitSystem>
    bool writeParticleData(const ParticleType* particles, S64 num_particles);
};

// CSV Output Writer
class CSVWriter : public OutputWriter {
private:
    std::ofstream file_;
    std::vector<std::string> column_names_;
    std::vector<std::string> column_units_;
    bool header_written_;

public:
    CSVWriter(const char* filename, const OutputMetadata& metadata, const CommInfo& comm_info)
        : OutputWriter(filename, metadata, comm_info), header_written_(false) {}

    ~CSVWriter() override {
        if (isOpen()) {
            close();
        }
    }

    bool open() override;
    void close() override;
    bool isOpen() const override { return file_.is_open(); }

    // Set column names and units
    void setColumnNames(const std::vector<std::string>& names) {
        column_names_ = names;
    }

    void setColumnUnits(const std::vector<std::string>& units) {
        column_units_ = units;
    }

    template<typename ParticleType, typename UnitSystem>
    bool writeParticleDataImpl(const ParticleType* particles, S64 num_particles);

protected:
    bool writeHeader() override;
    bool writeFooter() override;

private:
};

// HDF5 Output Writer (declaration only - implementation requires HDF5 library)
#ifdef PARTICLE_SIMULATOR_USE_HDF5
#include <hdf5.h>

class HDF5Writer : public OutputWriter {
private:
    hid_t file_id_;
    hid_t particle_group_id_;
    bool is_open_;

public:
    HDF5Writer(const char* filename, const OutputMetadata& metadata, const CommInfo& comm_info)
        : OutputWriter(filename, metadata, comm_info),
          file_id_(-1),
          particle_group_id_(-1),
          is_open_(false) {}

    ~HDF5Writer() override {
        if (isOpen()) {
            close();
        }
    }

    bool open() override;
    void close() override;
    bool isOpen() const override { return is_open_; }

    template<typename ParticleType, typename UnitSystem>
    bool writeParticleDataImpl(const ParticleType* particles, S64 num_particles);

protected:
    bool writeHeader() override;
    bool writeFooter() override;

private:
    bool writeAttribute(hid_t loc_id, const char* name, const char* value);
    bool writeAttribute(hid_t loc_id, const char* name, F64 value);
    bool writeAttribute(hid_t loc_id, const char* name, S64 value);
};
#endif

// Factory function to create appropriate writer
inline std::unique_ptr<OutputWriter> createOutputWriter(
    OutputFormat format,
    const char* filename,
    const OutputMetadata& metadata,
    const CommInfo& comm_info = CommInfo())
{
    switch (format) {
        case OutputFormat::CSV:
            return std::make_unique<CSVWriter>(filename, metadata, comm_info);
#ifdef PARTICLE_SIMULATOR_USE_HDF5
        case OutputFormat::HDF5:
            return std::make_unique<HDF5Writer>(filename, metadata, comm_info);
#endif
        default:
            PARTICLE_SIMULATOR_PRINT_ERROR("Unsupported output format");
            return nullptr;
    }
}

// Helper function to get file extension for format
inline const char* getFileExtension(OutputFormat format) noexcept {
    switch (format) {
        case OutputFormat::CSV:
            return ".csv";
        case OutputFormat::HDF5:
            return ".h5";
        case OutputFormat::Binary:
            return ".bin";
        case OutputFormat::ASCII:
            return ".txt";
        default:
            return "";
    }
}

// ===========================================
// Implementation of CSVWriter methods
// ===========================================

inline bool CSVWriter::open() {
    if (isOpen()) {
        return true;
    }

    // For MPI, append rank to filename
    std::string full_filename = filename_;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if (comm_info_.getNumberOfProc() > 1) {
        std::ostringstream oss;
        oss << filename_ << "_rank" << std::setfill('0') << std::setw(4) 
            << comm_info_.getRank() << ".csv";
        full_filename = oss.str();
    } else {
        full_filename += ".csv";
    }
#else
    full_filename += ".csv";
#endif

    file_.open(full_filename, std::ios::out);
    if (!file_.is_open()) {
        PARTICLE_SIMULATOR_PRINT_ERROR("Failed to open CSV file");
        std::cerr << "Filename: " << full_filename << std::endl;
        return false;
    }

    return true;
}

inline void CSVWriter::close() {
    if (isOpen()) {
        writeFooter();
        file_.close();
    }
}

inline bool CSVWriter::writeHeader() {
    if (!isOpen() || header_written_) {
        return false;
    }

    // Write metadata as comments
    file_ << "# Particle Simulation Output (CSV Format)" << std::endl;
    file_ << "# Time: " << std::scientific << std::setprecision(15) << metadata_.time << " " 
          << metadata_.time_unit << std::endl;
    file_ << "# Total Particles: " << metadata_.total_particles << std::endl;
    file_ << "# Unit System: " << metadata_.unit_system_name << std::endl;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    file_ << "# MPI Rank: " << metadata_.mpi_rank << " of " << metadata_.mpi_size << std::endl;
#endif
    if (!metadata_.description.empty()) {
        file_ << "# Description: " << metadata_.description << std::endl;
    }
    file_ << "#" << std::endl;

    // Write column headers
    if (!column_names_.empty()) {
        for (size_t i = 0; i < column_names_.size(); ++i) {
            if (i > 0) file_ << ",";
            file_ << column_names_[i];
        }
        file_ << std::endl;
    }

    // Write column units
    if (!column_units_.empty()) {
        file_ << "# Units: ";
        for (size_t i = 0; i < column_units_.size(); ++i) {
            if (i > 0) file_ << ",";
            file_ << column_units_[i];
        }
        file_ << std::endl;
    }

    header_written_ = true;
    return true;
}

inline bool CSVWriter::writeFooter() {
    if (!isOpen()) {
        return false;
    }
    file_ << "# End of file" << std::endl;
    return true;
}

// Particle data writer template specialization helper
template<typename ParticleType, typename UnitSystem>
inline bool CSVWriter::writeParticleDataImpl(const ParticleType* particles, S64 num_particles) {
    if (!isOpen()) {
        return false;
    }

    // Set up default column names if not already set
    if (column_names_.empty()) {
        column_names_ = {"id", "mass", "pos_x", "pos_y", "pos_z",
                        "vel_x", "vel_y", "vel_z", "dens", "eng", "pres"};
    }

    // Set up default unit strings
    if (column_units_.empty()) {
        using namespace UnitSystem;
        column_units_ = {
            "none",
            getUnitString<MassDim, UnitSystem>(),
            getUnitString<LengthDim, UnitSystem>(),
            getUnitString<LengthDim, UnitSystem>(),
            getUnitString<LengthDim, UnitSystem>(),
            getUnitString<VelocityDim, UnitSystem>(),
            getUnitString<VelocityDim, UnitSystem>(),
            getUnitString<VelocityDim, UnitSystem>(),
            getUnitString<DensityDim, UnitSystem>(),
            getUnitString<EnergyDim, UnitSystem>(),
            getUnitString<PressureDim, UnitSystem>()
        };
    }

    // Write header if not done yet
    if (!header_written_) {
        writeHeader();
    }

    // Write particle data with high precision
    file_ << std::scientific << std::setprecision(15);
    for (S64 i = 0; i < num_particles; ++i) {
        file_ << particles[i].id << ","
              << particles[i].mass << ","
              << particles[i].pos.x << ","
              << particles[i].pos.y << ","
              << particles[i].pos.z << ","
              << particles[i].vel.x << ","
              << particles[i].vel.y << ","
              << particles[i].vel.z << ","
              << particles[i].dens << ","
              << particles[i].eng << ","
              << particles[i].pres << std::endl;
    }

    return true;
}

// Template method pattern implementation for base class
template<typename ParticleType, typename UnitSystem>
inline bool OutputWriter::writeParticles(const ParticleType* particles, S64 num_particles) {
    if (!open()) {
        return false;
    }

    // Update metadata
    using namespace UnitSystem;
    metadata_.unit_system_name = getUnitSystemName<UnitSystem>();

    // Write data using derived class implementation
    bool success = writeParticleData<ParticleType, UnitSystem>(particles, num_particles);

    return success;
}

template<typename ParticleType, typename UnitSystem>
inline bool OutputWriter::writeParticleData(const ParticleType* particles, S64 num_particles) {
    // Dispatch to concrete implementation
    if (auto* csv_writer = dynamic_cast<CSVWriter*>(this)) {
        return csv_writer->writeParticleDataImpl<ParticleType, UnitSystem>(particles, num_particles);
    }
#ifdef PARTICLE_SIMULATOR_USE_HDF5
    if (auto* hdf5_writer = dynamic_cast<HDF5Writer*>(this)) {
        return hdf5_writer->writeParticleDataImpl<ParticleType, UnitSystem>(particles, num_particles);
    }
#endif
    return false;
}

#ifdef PARTICLE_SIMULATOR_USE_HDF5
// ===========================================
// Implementation of HDF5Writer methods
// ===========================================

inline bool HDF5Writer::open() {
    if (isOpen()) {
        return true;
    }

    // Create filename with rank for MPI
    std::string full_filename = filename_;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    if (comm_info_.getNumberOfProc() > 1) {
        std::ostringstream oss;
        oss << filename_ << "_rank" << std::setfill('0') << std::setw(4) 
            << comm_info_.getRank() << ".h5";
        full_filename = oss.str();
    } else {
        full_filename += ".h5";
    }
#else
    full_filename += ".h5";
#endif

    // Create HDF5 file
    file_id_ = H5Fcreate(full_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id_ < 0) {
        PARTICLE_SIMULATOR_PRINT_ERROR("Failed to create HDF5 file");
        std::cerr << "Filename: " << full_filename << std::endl;
        return false;
    }

    is_open_ = true;
    return true;
}

inline void HDF5Writer::close() {
    if (isOpen()) {
        writeFooter();
        if (particle_group_id_ >= 0) {
            H5Gclose(particle_group_id_);
            particle_group_id_ = -1;
        }
        if (file_id_ >= 0) {
            H5Fclose(file_id_);
            file_id_ = -1;
        }
        is_open_ = false;
    }
}

inline bool HDF5Writer::writeHeader() {
    if (!isOpen()) {
        return false;
    }

    // Create group for particles
    particle_group_id_ = H5Gcreate(file_id_, "/Particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (particle_group_id_ < 0) {
        PARTICLE_SIMULATOR_PRINT_ERROR("Failed to create particle group in HDF5");
        return false;
    }

    // Write metadata as attributes
    writeAttribute(file_id_, "Time", metadata_.time);
    writeAttribute(file_id_, "TotalParticles", metadata_.total_particles);
    writeAttribute(file_id_, "UnitSystem", metadata_.unit_system_name.c_str());
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    writeAttribute(file_id_, "MPIRank", static_cast<S64>(metadata_.mpi_rank));
    writeAttribute(file_id_, "MPISize", static_cast<S64>(metadata_.mpi_size));
#endif
    if (!metadata_.description.empty()) {
        writeAttribute(file_id_, "Description", metadata_.description.c_str());
    }

    return true;
}

inline bool HDF5Writer::writeFooter() {
    return true;
}

inline bool HDF5Writer::writeAttribute(hid_t loc_id, const char* name, const char* value) {
    hid_t attr_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(attr_type, strlen(value) + 1);
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate(loc_id, name, attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    
    herr_t status = H5Awrite(attr_id, attr_type, value);
    
    H5Aclose(attr_id);
    H5Sclose(attr_space);
    H5Tclose(attr_type);
    
    return status >= 0;
}

inline bool HDF5Writer::writeAttribute(hid_t loc_id, const char* name, F64 value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate(loc_id, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    
    herr_t status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
    
    H5Aclose(attr_id);
    H5Sclose(attr_space);
    
    return status >= 0;
}

inline bool HDF5Writer::writeAttribute(hid_t loc_id, const char* name, S64 value) {
    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate(loc_id, name, H5T_NATIVE_LLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    
    herr_t status = H5Awrite(attr_id, H5T_NATIVE_LLONG, &value);
    
    H5Aclose(attr_id);
    H5Sclose(attr_space);
    
    return status >= 0;
}

template<typename ParticleType, typename UnitSystem>
inline bool HDF5Writer::writeParticleDataImpl(const ParticleType* particles, S64 num_particles) {
    if (!isOpen()) {
        return false;
    }

    // Write header if not done yet
    if (particle_group_id_ < 0) {
        writeHeader();
    }

    using namespace UnitSystem;

    // Define dataset dimensions
    hsize_t dims[1] = {static_cast<hsize_t>(num_particles)};
    hid_t dataspace = H5Screate_simple(1, dims, NULL);

    // Helper lambda to write a dataset with units
    auto writeDataset = [&](const char* name, const char* unit_str, const F64* data) {
        hid_t dataset = H5Dcreate(particle_group_id_, name, H5T_NATIVE_DOUBLE, dataspace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset < 0) return false;

        herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        
        // Write unit as attribute
        writeAttribute(dataset, "Units", unit_str);
        
        H5Dclose(dataset);
        return status >= 0;
    };

    // Extract and write each field
    std::vector<S64> ids(num_particles);
    std::vector<F64> mass(num_particles), pos_x(num_particles), pos_y(num_particles), pos_z(num_particles);
    std::vector<F64> vel_x(num_particles), vel_y(num_particles), vel_z(num_particles);
    std::vector<F64> dens(num_particles), eng(num_particles), pres(num_particles);

    for (S64 i = 0; i < num_particles; ++i) {
        ids[i] = particles[i].id;
        mass[i] = particles[i].mass;
        pos_x[i] = particles[i].pos.x;
        pos_y[i] = particles[i].pos.y;
        pos_z[i] = particles[i].pos.z;
        vel_x[i] = particles[i].vel.x;
        vel_y[i] = particles[i].vel.y;
        vel_z[i] = particles[i].vel.z;
        dens[i] = particles[i].dens;
        eng[i] = particles[i].eng;
        pres[i] = particles[i].pres;
    }

    // Write ID dataset (integer)
    hid_t id_dataset = H5Dcreate(particle_group_id_, "ID", H5T_NATIVE_LLONG, dataspace,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(id_dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids.data());
    H5Dclose(id_dataset);

    // Write floating-point datasets with units
    writeDataset("Mass", getUnitString<MassDim, UnitSystem>(), mass.data());
    writeDataset("PositionX", getUnitString<LengthDim, UnitSystem>(), pos_x.data());
    writeDataset("PositionY", getUnitString<LengthDim, UnitSystem>(), pos_y.data());
    writeDataset("PositionZ", getUnitString<LengthDim, UnitSystem>(), pos_z.data());
    writeDataset("VelocityX", getUnitString<VelocityDim, UnitSystem>(), vel_x.data());
    writeDataset("VelocityY", getUnitString<VelocityDim, UnitSystem>(), vel_y.data());
    writeDataset("VelocityZ", getUnitString<VelocityDim, UnitSystem>(), vel_z.data());
    writeDataset("Density", getUnitString<DensityDim, UnitSystem>(), dens.data());
    writeDataset("Energy", getUnitString<EnergyDim, UnitSystem>(), eng.data());
    writeDataset("Pressure", getUnitString<PressureDim, UnitSystem>(), pres.data());

    H5Sclose(dataspace);

    return true;
}
#endif // PARTICLE_SIMULATOR_USE_HDF5

} // namespace Output
} // namespace ParticleSimulator
