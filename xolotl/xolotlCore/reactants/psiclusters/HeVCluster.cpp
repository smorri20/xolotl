// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHe), numV(numV) {

	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Update the composition map
	compositionMap[heType] = numHe;
	compositionMap[vType] = numV;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "HeV";

	// Compute the reaction radius
	reactionRadius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
			+ pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0) * numV)
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
			- pow(
					(3.0 * pow(xolotlCore::latticeConstant, 3.0))
							/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

}

HeVCluster::HeVCluster(const HeVCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;
}

HeVCluster::~HeVCluster() {
}

std::shared_ptr<Reactant> HeVCluster::clone() {
	std::shared_ptr<Reactant> reactant(new HeVCluster(*this));
	return reactant;
}

double HeVCluster::getGenByEm() {
	return 0;
}

double HeVCluster::getAnnByEm() {
	return 0;
}

void HeVCluster::createReactionConnectivity() {

	// Local Declarations
	auto psiNetwork = std::dynamic_pointer_cast<PSIClusterReactionNetwork>(
			network);
	auto props = psiNetwork->getProperties();
	int maxHeClusterSize = std::stoi(props["maxHeClusterSize"]);
	int maxHeVClusterSize = std::stoi(props["maxHeVClusterSize"]);
	std::vector<int> firstComposition, secondComposition;

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(getId());

	/* ----- (He_a)(V_b) + (He_c) --> [He_(a+c)]*(V_b) -----
	 * Helium absorption by HeV clusters that results
	 * in the production of this cluster.
	 */
	for (int z = 1; z <= maxHeClusterSize; z++) {
		// Get the first reactant
		firstComposition = psiNetwork->getCompositionVector(numHe - z, numV, 0);
		auto firstReactant = psiNetwork->getCompound(typeName, firstComposition);
		// Get the second reactant
		auto secondReactant = psiNetwork->get(heType, z);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant) {
			ClusterPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
			// Add the pair to the list
			reactingPairs.push_back(pair);
		}
	}

	/* ----- (He_a)(V_b) + V --> (He_a)[V_(b+1)] -----
	 * HeV clusters are also produced by single-vacancy absorption by another
	 * HeV cluster. In this case, (A*He)[(B-1)*V] produces the current cluster.
	 */
	firstComposition = psiNetwork->getCompositionVector(numHe, numV - 1, 0);
	auto firstReactant = psiNetwork->getCompound(typeName, firstComposition);
	auto secondReactant = psiNetwork->get(vType, 1);
	// Create a ReactingPair with the two reactants
	if (firstReactant && secondReactant) {
		ClusterPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);
		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	/* ----- (He_a) + (V_b) --> (He_a)(V_b) -----
	 * Helium-vacancy clustering that results
	 * in the production of this cluster.
	 */
	// Get the first reactant
	firstReactant = psiNetwork->get(heType, numHe);
	// Get the second reactant
	secondReactant = psiNetwork->get(vType, numV);
	// Create a ReactingPair with the two reactants
	if (firstReactant && secondReactant) {
		ClusterPair pair((PSICluster *) firstReactant, (PSICluster *) secondReactant);;
		// Add the pair to the list
		reactingPairs.push_back(pair);
	}

	/* ----- (He_a)*(V_b) + I_c  --> (He_a)*[V_(b-c)] -----
	 * Helium-vacancy clusters lose vacancies when they interact with
	 * interstitial clusters.
	 *
	 * We assume that the HeV and interstitial cluster can only
	 * interact if they would produce another HeV cluster, not single He.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the (outgoing) flux due to combination
	 * reactions.
	 */
	auto reactants = psiNetwork->getAll(iType);
	replaceInCompound(reactants, iType, vType);

	/* ---- (He_a)*(V_b) + He_c --> [He_(a+c)]*(V_b) ----
	 * HeV clusters can absorb helium clusters so long as they do not cross
	 * the max size limit.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	reactants = psiNetwork->getAll(heType);
	combineClusters(reactants, maxHeVClusterSize, typeName);

	/* ----- (He_a)*(V_b) + V --> (He_a)*[V_(b+1)] -----
	 * HeV clusters can absorb single vacancies.
	 *
	 * All of these clusters are added to the set of combining reactants
	 * because they contribute to the flux due to combination reactions.
	 */
	secondReactant = psiNetwork->get(vType, 1);
	if (secondReactant) {
		// Create a container for it
		std::vector<Reactant *> singleVInVector;
		singleVInVector.push_back(secondReactant);
		// Call the combination function even though there is only one cluster
		// because it handles all of the work to properly connect the three
		// clusters in the reaction.
		combineClusters(singleVInVector,maxHeVClusterSize,typeName);
	}

	// Set the references to the size one clusters
	heCluster = (PSICluster *) network->get(heType, 1);
	vCluster = (PSICluster *) network->get(vType, 1);
	iCluster = (PSICluster *) network->get(iType, 1);

	// Set the references for clusters one size smaller than this cluster,
	// starting with one less He.
	std::vector<int> compositionVec = { numHe - 1, numV, 0 };
	heVClusterLessHe = (PSICluster *) network->getCompound(heVType, compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, numV - 1, 0};
	heVClusterLessV = (PSICluster *) network->getCompound(heVType, compositionVec);

	return;
}

void HeVCluster::createDissociationConnectivity() {
	// Store the cluster with one less helium
	std::vector<int> compositionVec = { numHe - 1, numV, 0 };
	auto heVClusterLessHe = network->getCompound(typeName, compositionVec);
	// Store the cluster with one more helium
	compositionVec = { numHe + 1, numV, 0 };
	auto heVClusterMoreHe = network->getCompound(typeName, compositionVec);
	// Store the cluster with one less vacancy
	compositionVec = {numHe, numV - 1, 0};
	auto heVClusterLessV = network->getCompound(typeName, compositionVec);
	// Store the cluster with one more vacancy
	compositionVec = {numHe, numV + 1, 0};
	auto heVClusterMoreV = network->getCompound(typeName, compositionVec);

	// He Dissociation
	// (He_a)(V_b) --> [He_(a-1)](V_b) + He_1
	auto singleCluster = network->get(heType, 1);
	emitClusters(singleCluster, heVClusterLessHe);
	// [He_(a+1)](V_b) --> (He_a)(V_b) + He_1
	// Here it is important that heVClusterMoreHe is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreHe, singleCluster);

	// Vacancy Dissociation
	// (He_a)(V_b) --> He_(a)[V_(b-1)] + V_1
	singleCluster = network->get(vType, 1);
	emitClusters(singleCluster, heVClusterLessV);
	// He_(a)[V_(b+1)] --> (He_a)(V_b) + V_1
	// Here it is important that heVClusterMoreV is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterMoreV, singleCluster);

	// Trap mutation
	// (He_a)(V_b) --> He_(a)[V_(b+1)] + I_1
	singleCluster = network->get(iType, 1);
	emitClusters(singleCluster, heVClusterMoreV);
	// He_(a)[V_(b-1)] --> (He_a)(V_b) + I_1
	// Here it is important that heVClusterLessV is the first cluster
	// because it is the dissociating one.
	dissociateCluster(heVClusterLessV, singleCluster);

	return;
}

void HeVCluster::setTemperature(double temp) {

	// Call the base class version to set all of the basic quantities.
	PSICluster::setTemperature(temp);

	// Calculate the much easier f4 term... first
	f4 = calculateDissociationConstant(*this, *heCluster, temperature)
			+ calculateDissociationConstant(*this, *vCluster, temperature)
			+ calculateDissociationConstant(*this, *iCluster, temperature);

	return;
}
