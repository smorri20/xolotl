// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int _numHe, int _numV,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
    PSICluster(_network, registry, buildName(_numHe, _numV)),
        numHe(_numHe), numV(_numV) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Update the composition map
	composition[toCompIdx(Species::He)] = numHe;
	composition[toCompIdx(Species::V)] = numV;

	// Set the typename appropriately
	type = ReactantType::HeV;

	// Compute the reaction radius
	// It is the same formula for HeV clusters
	reactionRadius = xolotlCore::tungstenLatticeConstant
			* pow((3.0 * numV) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;

	return;
}

