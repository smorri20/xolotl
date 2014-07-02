// Includes
#include "VCluster.h"
#include <iostream>
#include <Constants.h>
#include <PSIClusterReactionNetwork.h>

using namespace xolotlCore;

VCluster::VCluster(int nV, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(nV, registry) {

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "V_" << size;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "V";

	// Update the composition map
	compositionMap["V"] = size;

	// Compute the reaction radius
	// FIXME Not right...
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant;

}

VCluster::~VCluster() {
}

std::shared_ptr<Reactant> VCluster::clone() {
	std::shared_ptr<Reactant> reactant(new VCluster(*this));
	return reactant;
}

void VCluster::createReactionConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species reaction
	PSICluster::createReactionConnectivity();

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Get all the He clusters from the network
	auto reactants = network->getAll(heType);
	// combineClusters handles He combining with V to form HeV
	combineClusters(reactants, heVType);

	// Single Vacancy absorption by HeV clusters
	// (He_a)(V_b) + V --> (He_a)[V_(b+1)]
	// Only if the size of this cluster is 1
	if (size == 1) {
		// Get all the HeV clusters from the network
		reactants = network->getAll(heVType);
		// combineClusters handles HeV combining with V to form HeV
		combineClusters(reactants, heVType);
	}

	// Vacancy-Interstitial annihilation
	// I_a + V_b
	//        --> I_(a-b), if a > b
	//        --> V_(b-a), if a < b
	//        --> 0, if a = b
	// Get all the I clusters from the network
	reactants = network->getAll(iType);
	// fillVWithI handles this reaction
	fillVWithI(iType, reactants);

	// Add the reactants to the reacting pairs for the cases where this cluster is
	// produced by the above reaction. This has to be checked for every interstitial.
	auto reactantsSize = reactants.size();
	for (int i = 0; i < reactantsSize; i++) {
		auto firstReactant = (PSICluster *) reactants[i];
		// Get the interstitial cluster that is bigger than the vacancy
		// and can form this cluster. V only results when it is bigger than I.
		auto secondReactant = (PSICluster *) network->get(typeName, firstReactant->getSize() + size);
		// Add to the reacting pairs if the second reactant exists
		if (secondReactant) {
			// Create the pair
			ClusterPair pair(firstReactant, secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	// Interstitial reduction by Vacancy absorption in HeI clusters
	// (He_a)*(I_b) + (V_c) --> (He_a)*[I_(b-c)]
	// Get all the HeI clusters from the network
	reactants = network->getAll(heIType);
	// replaceInCompound handles this reaction
	replaceInCompound(reactants, iType, vType);

	return;
}

void VCluster::createDissociationConnectivity() {
	// Call the function from the PSICluster class to take care of the single
	// species dissociation
	PSICluster::createDissociationConnectivity();

	// Specific case for the single species cluster
	if (size == 1) {
		// V dissociation of HeV cluster is handled here
		// (He_a)(V_b) --> (He_a)[V_(b-1)] + V_1
		// Get all the HeV clusters of the network
		auto allHeVReactants = network->getAll(heVType);
		for (int i = 0; i < allHeVReactants.size(); i++) {
			auto cluster = (PSICluster *) allHeVReactants.at(i);

			// (He_a)(V_b) is the dissociating one, (He_a)[V_(b-1)] is the one
			// that is also emitted during the dissociation
			auto comp = cluster->getComposition();
			std::vector<int> compositionVec = { comp[heType], comp[vType] - 1,
					comp[iType] };
			auto smallerReactant = network->getCompound(heVType, compositionVec);
			dissociateCluster(allHeVReactants.at(i), smallerReactant);
		}
	}

	return;
}
