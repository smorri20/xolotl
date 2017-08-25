#ifndef VCLUSTER_H
#define VCLUSTER_H

// Includes
#include "PSICluster.h"
#include <xolotlPerf.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of atomic vacancies.
 */
class VCluster: public PSICluster {

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    VCluster() = delete;

	/**
	 * The constructor. All VClusters must be initialized with a size.
	 *
	 * @param nV the number of atomic vacancies in the cluster
	 * @param registry The performance handler registry
	 */
	VCluster(int nV,
        IReactionNetwork& _network,
        std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

    /**
     * Copy constructor, deleted to prevent use.
     */
    VCluster(const VCluster& other) = delete;


	//! Destructor
	~VCluster() {
	}

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux() const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const;

};
//end class VCluster

} /* end namespace xolotlCore */

#endif
