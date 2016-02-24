#ifndef WSRIMADVECTIONHANDLER_H
#define WSRIMADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class WSRIMAdvectionHandler: public AdvectionHandler {

public:

	//! The Constructor
	WSRIMAdvectionHandler() {}

	//! The Destructor
	~WSRIMAdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for a
	 * (SRIM) tungsten material.
	 *
	 * @param network The network
	 */
	void initialize(PSIClusterReactionNetwork *network) {
		// Get all the reactants and their number
//		auto reactants = network->getAll();
//		int size = reactants->size();

		// Clear the index and sink strength vectors
		indexVector.clear();
		sinkStrengthVector.clear();

		return;
	}

};
//end class WSRIMAdvectionHandler

} /* end namespace xolotlCore */
#endif
