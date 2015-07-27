#ifndef ADVECTIONHANDLER_H
#define ADVECTIONHANDLER_H

// Includes
#include "IAdvectionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium clusters. It needs to have subclasses
 * that implement the initialize method.
 */
class AdvectionHandler: public IAdvectionHandler {
protected:

	//! The vector containing the indices of the advecting clusters
	std::vector<int> indexVector;

	//! The vector containing the value of the sink strength (called A) of the advecting clusters
	std::vector<double> sinkStrengthVector;

public:

	//! The Constructor
	AdvectionHandler() {}

	//! The Destructor
	~AdvectionHandler() {}

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	int getNumberOfAdvecting() {return indexVector.size();}

};
//end class AdvectionHandler

} /* end namespace xolotlCore */
#endif
