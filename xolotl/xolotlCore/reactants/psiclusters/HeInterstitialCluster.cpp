// Includes
#include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeInterstitialCluster::HeInterstitialCluster(int numHelium, int numInterstitial,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHelium), numI(numInterstitial) {

	// Set the cluster size as the sum of
	// the number of Helium and Interstitials
	size = numHe + numI;

	// Update the composition map
	compositionMap[heType] = numHe;
	compositionMap[iType] = numI;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "I_" << numI;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "HeI";

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0) * numI)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

}

HeInterstitialCluster::HeInterstitialCluster(const HeInterstitialCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numI = other.numI;
}

HeInterstitialCluster::~HeInterstitialCluster() {
}

std::shared_ptr<Reactant> HeInterstitialCluster::clone() {
	std::shared_ptr<Reactant> reactant = std::make_shared<HeInterstitialCluster>(*this);
	return reactant;
}

double HeInterstitialCluster::getGenByEm() {
	return 0;
}

double HeInterstitialCluster::getAnnByEm() {
	return 0;
}

void HeInterstitialCluster::createReactionConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	int maxVClusterSize = std::stoi(props["maxVClusterSize"]);
	int maxHeIClusterSize = std::stoi(props["maxHeIClusterSize"]);
	std::vector<int> firstComposition, secondComposition;

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	/* ----- (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)] -----
	 * This section adds the clusters that produce this cluster to the array
	 * for vacancy absorption by HeI.
	 */
	for (int c = 1; c <= maxVClusterSize; c++) { // hehe... c++!
		// Get the first reactant's composition and then retrieve it
		firstComposition = psiNetwork->getCompositionVector(numHe, 0, numI + c);
		auto firstReactant = psiNetwork->getCompound(typeName, firstComposition);
		// Set the second reactant
		auto secondReactant = psiNetwork->get(vType, c);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ClusterPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ----- (He_a)(I_b) + I --> (He_a)[I_(b+1)] -----
	 * This section adds the clusters that produce this cluster to the array
	 * for single species interstitial absorption by HeI.
	 *
	 * This section also handles the case of this cluster combining with a
	 * single interstitial to produce an HeI cluster with one addition
	 * interstitial.
	 */
	firstComposition = psiNetwork->getCompositionVector(numHe, 0, numI - 1);
	// Get those Reactants from the network
	auto firstReactant = psiNetwork->getCompound(typeName, firstComposition);
	auto secondReactant = psiNetwork->get(iType, 1);
	if (firstReactant && secondReactant) {
		// Create the Reacting Pair
		ClusterPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
		// Add single I to the list of clusters this one interacts with if it
		// doesn't violate the maximum size limit.
		setReactionConnectivity(secondReactant->getId());
		combiningReactants.push_back(secondReactant);
	}

	/* ----- (He_a)(I_b) + He_c --> [He_(a+c)](I_b)
	 * HeI clusters can absorb helium clusters of any size so long as the
	 * maximum size limit is not violated.
	 */
	auto reactants = psiNetwork->getAll(heType);
	combineClusters(reactants, maxHeIClusterSize, typeName);

	/* ----- (A*He)(B*I) + (C*V) --> (A*He)[(B-C)*I] -----
	 * This section adds the clusters that are produced by this cluster to the
	 * array for vacancy absorption by HeI.
	 */
	reactants = psiNetwork->getAll(vType);
	replaceInCompound(reactants, vType, iType);

	// Set the references to the size one clusters
	heCluster = (PSICluster *) network->get(heType, 1);
	vCluster = (PSICluster *) network->get(vType, 1);
	iCluster = (PSICluster *) network->get(iType, 1);

	// Set the references for clusters one size smaller than this cluster,
	// starting with one less He.
	std::vector<int> compositionVec = { numHe - 1, 0, numI };
	heIClusterLessHe = (PSICluster *) network->getCompound(heIType, compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, 0, numI - 1};
	heIClusterLessI = (PSICluster *) network->getCompound(heIType, compositionVec);

	return;
}

void HeInterstitialCluster::createDissociationConnectivity() {
	// Store the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, 0, numI };
	auto heIClusterLessHe = network->getCompound(typeName, compositionVec);
	// Store the cluster with one more helium
	compositionVec = { numHe + 1, 0, numI };
	auto heIClusterMoreHe = network->getCompound(typeName, compositionVec);
	// Store the cluster with one less interstitial
	compositionVec = {numHe, 0, numI - 1};
	auto heIClusterLessI = network->getCompound(typeName, compositionVec);
	// Store the cluster with one more interstitial
	compositionVec = {numHe, 0, numI + 1};
	auto heIClusterMoreI = network->getCompound(typeName, compositionVec);

	// He Dissociation
	// (He_a)(I_b) --> [He_(a-1)](I_b) + He_1
	auto singleCluster = network->get(heType, 1);
	emitClusters(singleCluster, heIClusterLessHe);
	// [He_(a+1)](V_b) --> (He_a)(V_b) + He_1
	// Here it is important that heVClusterMoreHe is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heIClusterMoreHe, singleCluster);

	// Interstitial Dissociation
	// (He_a)(I_b) --> He_(a)[I_(b-1)] + I_1
	singleCluster = network->get(iType, 1);
	emitClusters(singleCluster, heIClusterLessI);
	// He_(a)[I_(b+1)] --> (He_a)(I_b) + I_1
	// Here it is important that heIClusterMoreI is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heIClusterMoreI, singleCluster);

	return;
}

void HeInterstitialCluster::setTemperature(double temp) {

	// Call the base class version to set all of the basic quantities.
	PSICluster::setTemperature(temp);

	// Calculate the much easier f4 term... first
	f4 = calculateDissociationConstant(*this, *heCluster, temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);

	return;
}
