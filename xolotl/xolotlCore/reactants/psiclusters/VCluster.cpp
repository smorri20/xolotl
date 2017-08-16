// Includes
#include "VCluster.h"
#include <iostream>
#include <Constants.h>
#include <PSIClusterReactionNetwork.h>

using namespace xolotlCore;

VCluster::VCluster(int nV,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry) {
	// Set the size
	size = nV;
	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "V_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = vType;

	// Update the composition map
	compositionMap[vType] = size;

	// Compute the reaction radius
	// It is the same formula for HeV clusters
	reactionRadius = xolotlCore::tungstenLatticeConstant
			* pow((3.0 * size) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;

	return;
}

double VCluster::getEmissionFlux() const {
	// Initial declarations
	double flux = PSICluster::getEmissionFlux();

	// Compute the loss to dislocation sinks
	if (size < 5) {
		// k^2 * D * C
		flux += xolotlCore::sinkStrength * diffusionCoefficient * concentration;
	}

	return flux;
}

void VCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	PSICluster::getEmissionPartialDerivatives(partials);

	// Compute the loss to dislocation sinks
	if (size < 5) {
		// k^2 * D * C
		partials[id - 1] -= xolotlCore::sinkStrength * diffusionCoefficient;
	}

	return;
}
