#include "PSICluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

PSICluster::PSICluster() :
		Reactant() {
	// Set the reactant name appropriately
	name = "PSICluster";

	return;
}

PSICluster::PSICluster(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Reactant(registry) {
	// Set the reactant name appropriately
	name = "PSICluster";

	return;
}

// The copy constructor
PSICluster::PSICluster(PSICluster &other) :
		Reactant(other), reactingPairs(other.reactingPairs), combiningReactants(
				other.combiningReactants), dissociatingPairs(
				other.dissociatingPairs), emissionPairs(other.emissionPairs) {
	// Recompute all of the temperature-dependent quantities
	setTemperature(other.getTemperature());

	return;
}

void PSICluster::createProduction(std::shared_ptr<ProductionReaction> reaction,
		int a, int b, int c, int d) {
	// Create a cluster pair from the given reaction
	ClusterPair pair((PSICluster *) reaction->first,
			(PSICluster *) reaction->second);
	// Add it
	reactingPairs.push_back(pair);
	auto it = reactingPairs.end() - 1;

	// Add the reaction to the network
	reaction = network->addProductionReaction(reaction);
	// Link it to the pair
	(*it).reaction = reaction;

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
			secondVDistance = 0.0;
	if (reaction->first->getType() == PSISuperType) {
		auto super = (PSICluster *) reaction->first;
		firstHeDistance = super->getHeDistance(c);
		firstVDistance = super->getVDistance(d);
	}
	if (reaction->second->getType() == PSISuperType) {
		auto super = (PSICluster *) reaction->second;
		secondHeDistance = super->getHeDistance(c);
		secondVDistance = super->getVDistance(d);
	}
	(*it).a00 += 1.0;
	(*it).a10 += firstHeDistance;
	(*it).a20 += firstVDistance;
	(*it).a01 += secondHeDistance;
	(*it).a02 += secondVDistance;
	(*it).a11 += firstHeDistance * secondHeDistance;
	(*it).a12 += firstHeDistance * secondVDistance;
	(*it).a21 += firstVDistance * secondHeDistance;
	(*it).a22 += firstVDistance * secondVDistance;

	return;
}

void PSICluster::createCombination(std::shared_ptr<ProductionReaction> reaction,
		int a, int b) {
	// Look for the other cluster
	IReactant * secondCluster;
	if (reaction->first->getId() == id)
		secondCluster = reaction->second;
	else
		secondCluster = reaction->first;

	// Check if the reaction was already added
	std::vector<CombiningCluster>::reverse_iterator it;
	for (it = combiningReactants.rbegin(); it != combiningReactants.rend();
			it++) {
		if (secondCluster == (*it).combining) {
			break;
		}
	}
	if (it == combiningReactants.rend()) {
		// It was not already in so add it
		// Creates the combining cluster
		CombiningCluster combCluster((PSICluster *) secondCluster);
		// Add it
		combiningReactants.push_back(combCluster);
		it = combiningReactants.rbegin();

		// Create the corresponding production reaction
		auto newReaction = std::make_shared<ProductionReaction>((*it).combining,
				this);
		// Add it to the network
		newReaction = network->addProductionReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}

	// Update the coefficients
	double heDistance = 0.0, vDistance = 0.0;
	if (secondCluster->getType() == PSISuperType) {
		auto super = (PSICluster *) secondCluster;
		heDistance = super->getHeDistance(a);
		vDistance = super->getVDistance(b);
	}
	(*it).a0 += 1.0;
	(*it).a1 += heDistance;
	(*it).a2 += vDistance;

	return;
}

void PSICluster::createDissociation(
		std::shared_ptr<DissociationReaction> reaction, int a, int b, int c,
		int d) {
	// Look for the other cluster
	IReactant * emittedCluster;
	if (reaction->first->getId() == id)
		emittedCluster = reaction->second;
	else
		emittedCluster = reaction->first;

	// Check if the reaction was already added
	std::vector<ClusterPair>::reverse_iterator it;
	for (it = dissociatingPairs.rbegin(); it != dissociatingPairs.rend();
			it++) {
		if (reaction->dissociating == (*it).first
				&& emittedCluster == (*it).second) {
			break;
		}
	}
	if (it == dissociatingPairs.rend()) {
		// It was not already in so add it
		// Create the pair of them where it is important that the
		// dissociating cluster is the first one
		ClusterPair pair((PSICluster *) reaction->dissociating,
				(PSICluster *) emittedCluster);
		// Add it
		dissociatingPairs.push_back(pair);
		it = dissociatingPairs.rbegin();

		// Create the corresponding dissociation reaction
		auto newReaction = std::make_shared<DissociationReaction>((*it).first,
				(*it).second, this);
		// Add it to the network
		newReaction = network->addDissociationReaction(newReaction);
		// Link it to the pair
		(*it).reaction = newReaction;
	}

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction->dissociating->getType() == PSISuperType) {
		auto super = (PSICluster *) reaction->dissociating;
		firstHeDistance = super->getHeDistance(a);
		firstVDistance = super->getVDistance(b);
	}
	(*it).a00 += 1.0;
	(*it).a10 += firstHeDistance;
	(*it).a20 += firstVDistance;

	return;
}

void PSICluster::createEmission(std::shared_ptr<DissociationReaction> reaction,
		int a, int b, int c, int d) {
	// Create the pair of emitted clusters
	ClusterPair pair((PSICluster *) reaction->first,
			(PSICluster *) reaction->second);
	// Add the pair to the emission pair vector
	emissionPairs.push_back(pair);
	auto it = emissionPairs.end() - 1;

	// Add the reaction to the network
	reaction = network->addDissociationReaction(reaction);
	// Link it to the pair
	(*it).reaction = reaction;

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {
	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); ++it) {
		connectivity[*it - 1] = 1;
	}

	return connectivity;
}

std::vector<int> PSICluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network->getDOF());
}

std::vector<int> PSICluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network->getDOF());
}

void PSICluster::resetConnectivities() {
	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(heMomId);
	setDissociationConnectivity(heMomId);
	setReactionConnectivity(vMomId);
	setDissociationConnectivity(vMomId);

	// Loop on the effective reacting pairs
	for (auto it = reactingPairs.begin(); it != reactingPairs.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->id);
		setReactionConnectivity((*it).second->id);
		setReactionConnectivity((*it).first->heMomId);
		setReactionConnectivity((*it).second->heMomId);
		setReactionConnectivity((*it).first->vMomId);
		setReactionConnectivity((*it).second->vMomId);
	}

	// Loop on the effective combining reactants
	for (auto it = combiningReactants.begin(); it != combiningReactants.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).combining->id);
		setReactionConnectivity((*it).combining->heMomId);
		setReactionConnectivity((*it).combining->vMomId);
	}

	// Loop on the effective dissociating pairs
	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end();
			++it) {
		// The cluster is connecting to the dissociating cluster which
		// is the first one by definition
		setDissociationConnectivity((*it).first->id);
		setDissociationConnectivity((*it).first->heMomId);
		setDissociationConnectivity((*it).first->vMomId);
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	return;
}

void PSICluster::setReactionNetwork(
		const std::shared_ptr<IReactionNetwork> reactionNetwork) {
	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

double PSICluster::getDissociationFlux() const {
	// Initial declarations
	double flux = 0.0;
	PSICluster *dissociatingCluster = nullptr;

	// Loop over all dissociating clusters that form this cluster
	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end();
			++it) {
		// Get the dissociating clusters
		dissociatingCluster = (*it).first;
		double l0A = dissociatingCluster->getConcentration(0.0, 0.0);
		double lHeA = dissociatingCluster->getHeMomentum();
		double lVA = dissociatingCluster->getVMomentum();

		// Calculate the Dissociation flux
		flux += (*it).reaction->kConstant
				* ((*it).a00 * l0A + (*it).a10 * lHeA + (*it).a20 * lVA);
	}

	// Return the flux
	return flux;
}

double PSICluster::getEmissionFlux() const {
	// Initial declarations
	double flux = 0.0;

	// Loop over all the pairs
	for (auto it = emissionPairs.begin(); it != emissionPairs.end(); ++it) {
		// Update the flux
		flux += (*it).reaction->kConstant;
	}

	return flux * concentration;
}

double PSICluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Loop over all the reacting pairs
	for (auto it = reactingPairs.begin(); it != reactingPairs.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration(0.0, 0.0);
		double l0B = secondReactant->getConcentration(0.0, 0.0);
		double lHeA = firstReactant->getHeMomentum();
		double lHeB = secondReactant->getHeMomentum();
		double lVA = firstReactant->getVMomentum();
		double lVB = secondReactant->getVMomentum();
		// Update the flux
		flux += (*it).reaction->kConstant
				* ((*it).a00 * l0A * l0B + (*it).a01 * l0A * lHeB
						+ (*it).a02 * l0A * lVB + (*it).a10 * lHeA * l0B
						+ (*it).a11 * lHeA * lHeB + (*it).a12 * lHeA * lVB
						+ (*it).a20 * lVA * l0B + (*it).a21 * lVA * lHeB
						+ (*it).a22 * lVA * lVB);
	}

	// Return the production flux
	return flux;
}

double PSICluster::getCombinationFlux() const {
	// Local declarations
	double flux = 0.0;
	PSICluster *combiningCluster = nullptr;

	// Loop over all possible clusters
	for (auto it = combiningReactants.begin(); it != combiningReactants.end();
			++it) {
		// Get the cluster that combines with this one
		combiningCluster = (*it).combining;
		double l0B = combiningCluster->getConcentration(0.0, 0.0);
		double lHeB = combiningCluster->getHeMomentum();
		double lVB = combiningCluster->getVMomentum();
		// Calculate the combination flux
		flux += (*it).reaction->kConstant
				* ((*it).a0 * l0B + (*it).a1 * lHeB + (*it).a2 * lVB);
	}

	return flux * concentration;
}

std::vector<double> PSICluster::getPartialDerivatives() const {
	// Local Declarations
	std::vector<double> partials(network->getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return partials;
}

void PSICluster::getPartialDerivatives(std::vector<double> & partials) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSICluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	double value = 0.0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	for (auto it = reactingPairs.begin(); it != reactingPairs.end(); ++it) {
		// Get the two reacting clusters
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		double l0A = firstReactant->getConcentration(0.0, 0.0);
		double l0B = secondReactant->getConcentration(0.0, 0.0);
		double lHeA = firstReactant->getHeMomentum();
		double lHeB = secondReactant->getHeMomentum();
		double lVA = firstReactant->getVMomentum();
		double lVB = secondReactant->getVMomentum();

		// Compute the contribution from the first part of the reacting pair
		value = (*it).reaction->kConstant;
		index = firstReactant->id - 1;
		partials[index] += value
				* ((*it).a00 * l0B + (*it).a01 * lHeB + (*it).a02 * lVB);
		index = firstReactant->heMomId - 1;
		partials[index] += value
				* ((*it).a10 * l0B + (*it).a11 * lHeB + (*it).a12 * lVB);
		index = firstReactant->vMomId - 1;
		partials[index] += value
				* ((*it).a20 * l0B + (*it).a21 * lHeB + (*it).a22 * lVB);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant->id - 1;
		partials[index] += value
				* ((*it).a00 * l0A + (*it).a10 * lHeA + (*it).a20 * lVA);
		index = secondReactant->heMomId - 1;
		partials[index] += value
				* ((*it).a01 * l0A + (*it).a11 * lHeA + (*it).a21 * lVA);
		index = secondReactant->vMomId - 1;
		partials[index] += value
				* ((*it).a02 * l0A + (*it).a12 * lHeA + (*it).a22 * lVA);
	}

	return;
}

void PSICluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int otherIndex = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	for (auto it = combiningReactants.begin(); it != combiningReactants.end();
			++it) {
		cluster = (PSICluster *) (*it).combining;
		double l0B = cluster->getConcentration(0.0, 0.0);
		double lHeB = cluster->getHeMomentum();
		double lVB = cluster->getVMomentum();

		// Remember that the flux due to combinations is OUTGOING (-=)!
		// Compute the contribution from this cluster
		partials[id - 1] -= (*it).reaction->kConstant
				* ((*it).a0 * l0B + (*it).a1 * lHeB + (*it).a2 * lVB);
		// Compute the contribution from the combining cluster
		value = (*it).reaction->kConstant * concentration;
		otherIndex = cluster->id - 1;
		partials[otherIndex] -= value * (*it).a0;
		otherIndex = cluster->heMomId - 1;
		partials[otherIndex] -= value * (*it).a1;
		otherIndex = cluster->vMomId - 1;
		partials[otherIndex] -= value * (*it).a2;
	}

	return;
}

void PSICluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	for (auto it = dissociatingPairs.begin(); it != dissociatingPairs.end();
			++it) {
		// Get the dissociating cluster
		cluster = (*it).first;
		value = (*it).reaction->kConstant;
		index = cluster->id - 1;
		partials[index] += value * (*it).a00;
		index = cluster->heMomId - 1;
		partials[index] += value * (*it).a10;
		index = cluster->vMomId - 1;
		partials[index] += value * (*it).a20;
	}

	return;
}

void PSICluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	for (auto it = emissionPairs.begin(); it != emissionPairs.end(); ++it) {
		// Modify the partial derivative. Remember that the flux
		// due to emission is OUTGOING (-=)!
		index = id - 1;
		partials[index] -= (*it).reaction->kConstant;
	}

	return;
}

void PSICluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void PSICluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double PSICluster::getLeftSideRate() const {
	// Initialize the rate and the cluster pointer
	double totalRate = 0.0;
	PSICluster *cluster = nullptr;

	// Loop on the combining reactants
	for (int i = 0; i < combiningReactants.size(); i++) {
		cluster = (PSICluster *) combiningReactants[i].combining;
		// Add the rate to the total rate
		totalRate += combiningReactants[i].reaction->kConstant
				* cluster->concentration;
	}

	// Loop on the emission pairs
	for (int i = 0; i < emissionPairs.size(); i++) {
		// Add the rate to the total rate
		totalRate += emissionPairs[i].reaction->kConstant;
	}

	return totalRate;
}

std::vector<int> PSICluster::getConnectivity() const {
	int connectivityLength = network->getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector is an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnVector[i] || dissociationConnVector[i];
	}

	return connectivity;
}
