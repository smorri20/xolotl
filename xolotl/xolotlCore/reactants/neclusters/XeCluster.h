#ifndef XECLUSTER_H
#define XECLUSTER_H

// Includes
#include "NECluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of helium.
 */
class XeCluster: public NECluster {

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    XeCluster() = delete;

	/**
	 * The constructor. All XeClusters must be initialized with a size.
	 *
	 * @param nXe the number of helium atoms in the cluster
	 * @param registry The performance handler registry
	 */
	XeCluster(int nXe,
        IReactionNetwork& _network,
        std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

    /**
     * Copy constructor, deleted to prevent use.
     */
    XeCluster(const XeCluster& other) = delete;

	/**
	 * Destructor
	 */
	~XeCluster() {}

}; //end class XeCluster

} /* namespace xolotlCore */
#endif
