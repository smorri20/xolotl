#ifndef HEVCLUSTER_H
#define HEVCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <map>

namespace xolotlCore {

/**
 *  A cluster composed of helium and vacancies
 */
class HeVCluster: public PSICluster {

private:
    // TODO do we need to keep these species counts here, 
    // since they are in the composition?

	//! The number of helium atoms in this cluster.
	int numHe;

	//! The number of atomic vacancies in this cluster.
	int numV;


	static
    std::string buildName(IReactant::SizeType nHe, IReactant::SizeType nV) {
        std::stringstream nameStream;
        nameStream << "He_" << nHe << "V_" << nV;
        return nameStream.str();
    }

public:

    /**
     * Default constructor, deleted because we require info to construct.
     */
    HeVCluster() = delete;

	/**
	 * The constructor. All HeVClusters must be initialized with a map
	 * that describes the species of which the cluster is composed. The map
	 * should contain as its keys the names of the species and the sizes of the
	 * species as its values. The names of the species must be one of
	 * {He,V}.
	 *
	 * @param numHe The number of helium atoms in this cluster
	 * @param numV The number of vacancies in this cluster
     * @param _network The network the cluster will belong to.
	 * @param registry The performance handler registry
	 */
	HeVCluster(int numHe, int numV,
            IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	HeVCluster(HeVCluster &other) = delete;

	//! Destructor
	~HeVCluster() {}

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {return true;}

};
//end class HeVCluster

} /* end namespace xolotlCore */
#endif
