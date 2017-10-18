#ifndef HCLUSTER_H
#define HCLUSTER_H

// Includes
#include "PSICluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of hydrogen.
class HCluster: public PSICluster {

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    HCluster() = delete;

	/**
	 * The constructor. All HClusters must be initialized with a size.
	 * @param nH the number of hydrogen atoms in the cluster
	 */
	HCluster(int nH,
            IReactionNetwork& _network,
            std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
            const std::string& _name = "Hydrogen") :
        PSICluster(_network, registry, _name) {

        // Set the size appropriately
        size = nH;
    }

    /**
     * Copy constructor, deleted to prevent use.
     */
    HCluster(const HCluster& other) = delete;

	//! The destructor
	~HCluster() {
    }

};
//end class HCluster

} /* end namespace xolotlCore */
#endif
