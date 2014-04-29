#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

namespace HDF5Utils {

	/**
	 * Create the HDF5 file with the needed structure.
	 */
	void initializeFile(int timeStep, int networkSize);

	/**
	 * Fill the header.
	 */
	void fillHeader(int physicalDim, int refinement, double time, double deltaTime);

	/**
	 * Fill the network dataset.
	 */
	void fillNetwork(std::shared_ptr<PSIClusterReactionNetwork> network);

	/**
	 * Fill the concentration dataset at a specific grid point.
	 */
	void fillConcentrations(double * concArray, double position);

	/**
	 * Add the data to the file and close it.
	 */
	void finalizeFile();

};

} /* namespace xolotlCore */
#endif
