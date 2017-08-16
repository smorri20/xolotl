#ifndef DCLUSTER_H
#define DCLUSTER_H

// Includes
#include "HCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of deuterium.
class DCluster : public HCluster {

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    DCluster() = delete;

	//! The constructor. All DClusters must be initialized with a size.
	DCluster(int nD,
            IReactionNetwork& _network,
            std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! Destructor
	~DCluster();

};
//end class DCluster

} /* end namespace xolotlCore */
#endif
