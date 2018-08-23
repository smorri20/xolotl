#ifndef FRANKCLUSTER_H
#define FRANKCLUSTER_H

// Includes
#include "AlloyCluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial.
 */
class FrankCluster: public AlloyCluster {

private:

	/**
	 * The default constructor is private because NEClusters must always be
	 * initialized with a size and performance handler registry
	 */
	FrankCluster() :
		AlloyCluster() {}

public:

	/**
	 * The constructor. All FrankClusters must be initialized with a size.
	 *
	 * @param n The size of the cluster
	 * @param registry The performance handler registry
	 */
	FrankCluster(int n, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Destructor
	 */
	~FrankCluster() {}

	/**
	 * Returns a reactant created using the copy constructor
	 */
	virtual std::shared_ptr<IReactant> clone() {
		return std::shared_ptr<IReactant> (new FrankCluster(*this));
	}

}; //end class FrankCluster

} /* namespace xolotlCore */
#endif
