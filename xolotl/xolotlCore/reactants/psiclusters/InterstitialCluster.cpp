#include "InterstitialCluster.h"
#include <Constants.h>
#include <iostream>

using namespace xolotlCore;

InterstitialCluster::InterstitialCluster(int nI,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry) {
	// Set the size
	size = nI;
	// Update the composition map
	compositionMap[iType] = size;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "I_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = iType;

	// Compute the reaction radius
	double EightPi = 8.0 * xolotlCore::pi;
	reactionRadius = xolotlCore::tungstenLatticeConstant
			* pow((3.0 / EightPi) * size, (1.0 / 3.0));

	return;
}

double InterstitialCluster::getEmissionFlux() const {
	// Initial declarations
	double flux = PSICluster::getEmissionFlux();

//	// Compute the loss to dislocation sinks
//	if (size < 2) {
//		// bias * k^2 * D * C
//		flux += sinkBias * sinkStrength * diffusionCoefficient * concentration;
//	}

	return flux;
}

void InterstitialCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	PSICluster::getEmissionPartialDerivatives(partials);

//	// Compute the loss to dislocation sinks
//	if (size < 2) {
//		// bias * k^2 * D * C
//		partials[id - 1] -= sinkBias * sinkStrength * diffusionCoefficient;
//	}

	return;
}
