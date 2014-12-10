#ifndef DUMMYDIFFUSIONHANDLER_H
#define DUMMYDIFFUSIONHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile cluster. Here it is a dummy class,
 * meaning that it should not do anything.
 */
class DummyDiffusionHandler: public DiffusionHandler {
public:

	//! The Constructor
	DummyDiffusionHandler() {}

	//! The Destructor
	~DummyDiffusionHandler() {}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped it
	 * won't be possible to set the partials for the diffusion.
	 *
	 * We don't want any cluster to diffuse, so nothing is set to 1 in ofill, and no index
	 * is added to the vector.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the diffusing clusters
	 */
	void initializeOFill(PSIClusterReactionNetwork *network, int *ofill) {
		// Clear the index vector
		indexVector.clear();

		// And don't do anything else
		return;
	}

};
//end class DummyDiffusionHandler

} /* end namespace xolotlCore */
#endif
