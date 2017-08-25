#ifndef TCLUSTER_H
#define TCLUSTER_H

// Includes
#include "HCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of tritium.
class TCluster : public HCluster {

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    TCluster() = delete;

	//! The constructor. All TClusters must be initialized with a size.
	TCluster(int nT,
            IReactionNetwork& _network,
            std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

    /**
     * Copy constructor, deleted to prevent use.
     */
    TCluster(const TCluster& other) = delete;

	//! Destructor
	~TCluster();

};
//end class TCluster

} /* end namespace xolotlCore */
#endif
