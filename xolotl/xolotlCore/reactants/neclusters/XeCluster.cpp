// Includes
#include "XeCluster.h"
#include "NEClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

XeCluster::XeCluster(int nXe,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		NECluster(_network, registry, buildName(nXe)) {

	// Set the size
	size = nXe;
	// Update the composition map
	composition[toCompIdx(Species::Xe)] = size;

	// Set the typename appropriately
	type = ReactantType::Xe;

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;
	reactionRadius = 1.05 * pow((3.0 * 85.0 * (double) size) / FourPi, (1.0/3.0)) / 10.0;
	if (size == 1) reactionRadius = 0.3;

	return;
}
