#ifndef HECLUSTER_H
#define HECLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class HeCluster: public PSICluster {

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    HeCluster() = delete;

	/**
	 * The constructor. All HeClusters must be initialized with a size.
	 *
	 * @param nHe the number of helium atoms in the cluster
	 * @param registry The performance handler registry
	 */
	HeCluster(int nHe,
            IReactionNetwork& _network,
            std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

    /**
     * Copy constructor, deleted to prevent use.
     */
    HeCluster(const HeCluster& other) = delete;

	/**
	 * Destructor
	 */
	~HeCluster() {}

}; //end class HeCluster

} /* namespace xolotlCore */
#endif
