// Includes
#include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(int numHelium, int numInterstitial,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry), numHe(numHelium), numI(numInterstitial) {
	// Set the cluster size as the sum of
	// the number of Helium and Interstitials
	size = numHe + numI;

	// Update the composition map
	composition[toCompIdx(Species::He)] = numHe;
	composition[toCompIdx(Species::I)] = numI;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "I_" << numI;
	name = nameStream.str();
	// Set the typename appropriately
	type = ReactantType::HeI;

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0) * numI)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::tungstenLatticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

	return;
}


