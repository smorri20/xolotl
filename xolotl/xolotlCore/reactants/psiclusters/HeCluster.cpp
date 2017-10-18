// Includes
#include "HeCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

HeCluster::HeCluster(int nHe, IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(nHe)) {

	// Set the size
	size = nHe;
	// Update the composition map
	composition[toCompIdx(Species::He)] = size;

	// Set the typename appropriately
	type = ReactantType::He;

	// Compute the reaction radius
	double FourPi = 4.0 * xolotlCore::pi;
	double aCubed = pow(xolotlCore::tungstenLatticeConstant, 3);
	double termOne = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * size,
			(1.0 / 3.0));
	double termTwo = pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
	reactionRadius = 0.3 + termOne - termTwo;

	// Bounds on He and V
	heBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>(size),
			static_cast<IReactant::SizeType>(size+1));
	vBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>(0),
			static_cast<IReactant::SizeType>(1));

	return;
}
