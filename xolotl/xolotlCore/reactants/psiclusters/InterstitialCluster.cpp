#include "InterstitialCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(nI)) {

	// Set the size
	size = nI;
	// Update the composition map
	composition[toCompIdx(Species::I)] = size;

	// Set the typename appropriately
	type = ReactantType::I;

	// Compute the reaction radius
	double EightPi = 8.0 * xolotlCore::pi;
	reactionRadius = xolotlCore::tungstenLatticeConstant
			* pow((3.0 / EightPi) * size, (1.0 / 3.0));

	return;
}

double InterstitialCluster::getEmissionFlux() const {
	// Initial declarations
	double flux = PSICluster::getEmissionFlux();

	// Compute the loss to dislocation sinks
	if (size < 2) {
		// bias * k^2 * D * C
		flux += sinkBias * sinkStrength * diffusionCoefficient * concentration;
	}

	return flux;
}

void InterstitialCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	PSICluster::getEmissionPartialDerivatives(partials);

	// Compute the loss to dislocation sinks
	if (size < 2) {
		// bias * k^2 * D * C
		partials[id - 1] -= sinkBias * sinkStrength * diffusionCoefficient;
	}

	return;
}
