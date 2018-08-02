#include "xolotlCore/io/TridynFile.h"

namespace xolotlCore {

const std::string TridynFile::lastTimestepAttrName = "lastTimeStep";

TridynFile::TridynFile(fs::path path, AccessMode mode, MPI_Comm comm)
  : HDF5File(path, mode, comm, true) {

    // Ensure we have a last timestep attribute.
    ScalarDataSpace lastTimestepAttrDSpace;
    Attribute<int> lastTimestepAttr(*this,
                                    lastTimestepAttrName,
                                    lastTimestepAttrDSpace);

    std::cout << "created TridynFile, id: " << getId() << std::endl;
}


// Version for 1D problems.
// Our data is thus 2D (one dim for gridpoints, one dim for data values).
void TridynFile::writeTimestep(int timestep,
                                uint32_t numGridpoints,
                                uint32_t offset,
                                const ConcsType1D& data) const {

    // Define the overall dataset.
    SimpleDataSpace<2>::Dimensions dsetDims =
                        { (hsize_t)numGridpoints, numValsPerGridpoint };
    SimpleDataSpace<2> dsetSpace(dsetDims);
    std::ostringstream dsetNameStr;
    dsetNameStr << "concs_" << timestep;
    DataSet<ValueType> dset(*this, dsetNameStr.str(), dsetSpace);

	// Write the dataset.
    dset.parWrite2D<numValsPerGridpoint>(getComm(), offset, data);

    // Update our idea of the last timestep written.
    Attribute<int> lastTimestepAttr(*this, lastTimestepAttrName);
    lastTimestepAttr.setTo(timestep);
}

} // namespace xolotlCore

