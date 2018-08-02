#ifndef XCORE_TRIDYNFILE_H
#define XCORE_TRIDYNFILE_H

#include <string>
#include <vector>
#include <tuple>
#include "xolotlCore/io/HDF5File.h"
#include "xolotlCore/io/HDF5Exception.h"


namespace xolotlCore {


// Class for reading and writing a TRIDYN checkpoint file.
class TridynFile : public HDF5File {
private:
    using ValueType = double;
    using TimestepType = int;

    //! Name of the last timestep attribute.
    static const std::string lastTimestepAttrName;

public:
	static constexpr auto numSpecies = 5;	// He, D, T, V, I
	static constexpr auto numValsPerGridpoint = numSpecies + 1;	// for the grid
	using ConcsType1D = DataSet<ValueType>::DataType2D<numValsPerGridpoint>;

    /**
     * Create/initialize or open a TRIDYN checkpoint file.
     *
     * @param path Path of file to create or open.
     * @param mode Access mode for file.
     * @param comm The MPI communicator used to access the file.
     */
    TridynFile(fs::path path,
            AccessMode mode = AccessMode::CreateOrTruncateIfExists,
            MPI_Comm comm = MPI_COMM_WORLD);


    /**
     * Write concentrations associated with a timestep when running
     * a 1D problem to our file.
     *
     * @param timestep The time step to associate with the concentrations.
     * @param numGridpointsToWrite Total number of gridpoints in the 1D problem.
     * @param myOffset Our offset within the overall space.
     * @param concs The concentration data for the gridpoints that my 
     *              process owns.
     */
    void writeTimestep(TimestepType timestep,
                        uint32_t numGridpointsToWrite,
                        uint32_t myOffset,
                        const ConcsType1D& concs) const;
};

} /* namespace xolotlCore */

#endif // XCORE_TRIDYNFILE_H

