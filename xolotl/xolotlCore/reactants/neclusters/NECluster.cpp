#include "NECluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <iostream>
#include <limits>
#include <MathUtils.h>

using namespace xolotlCore;

NECluster::NECluster() :
		Reactant() {
	// Set the reactant name appropriately
	name = "NECluster";

	return;
}

NECluster::NECluster(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Reactant(registry) {
	// Set the reactant name appropriately
	name = "NECluster";

	return;
}

// The copy constructor
NECluster::NECluster(NECluster &other) :
		Reactant(other),
		reactingPairs(other.reactingPairs),
		combiningReactants(other.combiningReactants),
		dissociatingPairs(other.dissociatingPairs),
		emissionPairs(other.emissionPairs) {
	// Recompute all of the temperature-dependent quantities
	setTemperature(other.getTemperature());

	return;
}

void NECluster::createReactionConnectivity() {
	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);

	// This cluster is always X_a

	// Initial declarations
	int firstSize = 0, secondSize = 0;

	// Single species clustering producing this cluster
	// X_(a-i) + X_i --> X_a
	for (firstSize = 1; firstSize <= 1; firstSize++) {
		// Set the size of the second reactant
		secondSize = size - firstSize;
		// Get the first and second reactants for the reaction
		auto firstReactant = (NECluster *) network->get(typeName, firstSize);
		auto secondReactant = (NECluster *) network->get(typeName, secondSize);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant &&
				(firstReactant->diffusionFactor > 0.0 || secondReactant->diffusionFactor > 0.0)) {
			// The reaction constant will be computed later, it is set to 0.0 for now
			ClusterPair pair(firstReactant, secondReactant, 0.0);
			// Add the pair to the list
			reactingPairs.push_back(pair);
			// Setup the connectivity array
			setReactionConnectivity(firstReactant->id);
			setReactionConnectivity(secondReactant->id);
		}
	}

	// Single species clustering
	// X_a + X_b --> X_(a+b)
	if (size == 1) {
		auto reactants = network->getAll(typeName);
		// combineClusters handles everything for this type of reaction
		NECluster::combineClusters(reactants, typeName);
	}
	else {
		std::vector<IReactant *> reactants;
		reactants.push_back(network->get(typeName, 1));
		// combineClusters handles everything for this type of reaction
		NECluster::combineClusters(reactants, typeName);
	}

	return;
}

void NECluster::createDissociationConnectivity() {
	// This cluster is always X_a

	// Single species dissociation
	// X_a --> X_(a-1) + X
	auto smallerReactant = (NECluster *) network->get(typeName, size - 1);
	auto singleCluster = (NECluster *) network->get(typeName, 1);
	emitClusters(singleCluster, smallerReactant);

	// Single species dissociation producing this cluster
	// X_(a+1) --> X_a + X
	auto biggerReactant = (NECluster *) network->get(typeName, size + 1);
	dissociateCluster(biggerReactant, singleCluster);

	// Specific case for the single size cluster
	// for a = 1
	if (size == 1) {
		// all the cluster of the same type dissociate into it
		auto allSameTypeReactants = network->getAll(typeName);
		for (unsigned int i = 0; i < allSameTypeReactants.size(); i++) {
			auto cluster = (NECluster *) allSameTypeReactants[i];
			// X_1 cannot dissociate and X_2 --> X + X was already
			// counted in the previous step
			if (cluster->size < 3)
				continue;
			// X_b is the dissociating one, X_(b-a) is the one
			// that is also emitted during the dissociation
			smallerReactant = (NECluster *) network->get(typeName,
					cluster->size - 1);
			dissociateCluster(cluster, smallerReactant);
		}
	}

	return;
}

double NECluster::calculateReactionRateConstant(
		const NECluster & firstReactant, const NECluster & secondReactant) const {
	// Get the reaction radii
	double r_first = firstReactant.reactionRadius;
	double r_second = secondReactant.reactionRadius;

	// Get the diffusion coefficients
	double firstDiffusion = firstReactant.diffusionCoefficient;
	double secondDiffusion = secondReactant.diffusionCoefficient;

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi * (r_first + r_second)
					* (firstDiffusion + secondDiffusion);
	return k_plus;
}

double NECluster::calculateDissociationConstant(
		const NECluster & dissociatingCluster,
		const NECluster & singleCluster, const NECluster & secondCluster) const {
	//! NEEDS TO BE UPDATED
	double atomicVolume = 0.5 * xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant * xolotlCore::uraniumDioxydeLatticeConstant;

	// Calculate the Reaction Rate Constant
	double kPlus = calculateReactionRateConstant(singleCluster, secondCluster);

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(dissociatingCluster, singleCluster,
			secondCluster);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

double NECluster::computeBindingEnergy(const NECluster & dissociatingCluster,
		const NECluster & singleCluster,
		const NECluster & secondCluster) const {
	// for the dissociation A --> B + C we need A binding energy
	// E_b(A) = E_f(B) + E_f(C) - E_f(A) where E_f is the formation energy
	double bindingEnergy = singleCluster.formationEnergy
			+ secondCluster.formationEnergy
			- dissociatingCluster.formationEnergy;
	return bindingEnergy;
}

void NECluster::dissociateCluster(NECluster * dissociatingCluster,
		NECluster * emittedCluster) {
	// Test if the dissociatingCluster and the emittedCluster exist
	if (dissociatingCluster && emittedCluster) {
		// Create the pair of them where it is important that the
		// dissociating cluster is the first one
		// The dissociating constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(dissociatingCluster, emittedCluster, 0.0);
		// Add the pair to the dissociating pair vector
		// The connectivity is handled in emitCluster
		dissociatingPairs.push_back(pair);

		// Take care of the connectivity
		setDissociationConnectivity(dissociatingCluster->id);

		// Add it to the list again if it the same as the other emitted cluster
		if (id == emittedCluster->id) dissociatingPairs.push_back(pair);
	}

	return;
}

void NECluster::emitClusters(NECluster * firstEmittedCluster,
		NECluster * secondEmittedCluster) {
	// Test if the emitted clusters exist
	if (firstEmittedCluster && secondEmittedCluster) {
		// Connect this cluster to itself since any reaction will affect it
		setDissociationConnectivity(id);

		// Add the pair of emitted clusters to the vector of emissionPairs
		// The first cluster is the size one one
		// The dissociating constant will be computed later, it is set to 0.0 for now
		ClusterPair pair(firstEmittedCluster, secondEmittedCluster, 0.0);
		emissionPairs.push_back(pair);
	}

	return;
}

void NECluster::combineClusters(std::vector<IReactant *> & reactants,
		const std::string& productName) {
	// Initial declarations
	std::map<std::string, int> myComposition = getComposition(),
			secondComposition;
	int numXe, numV, numI, secondNumXe, secondNumV, secondNumI, productSize;
	std::vector<int> compositionSizes { 0, 0, 0 };
	NECluster *productCluster, *secondCluster;
	// Setup the composition variables for this cluster
	numXe = myComposition[xeType];
	numV = myComposition[vType];
	numI = myComposition[iType];

	int reactantVecSize = reactants.size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second reactant, its composition and its index
		secondCluster = (NECluster *) reactants[i];
		secondComposition = secondCluster->getComposition();
		secondNumXe = secondComposition[xeType];
		secondNumV = secondComposition[vType];
		secondNumI = secondComposition[iType];
		// Compute the product size
		productSize = size + secondCluster->size;
		// Get and handle product for compounds
		if (productName == heVType || productName == heIType) {
			// Modify the composition vector
			compositionSizes[0] = numXe + secondNumXe;
			compositionSizes[1] = numV + secondNumV;
			compositionSizes[2] = numI + secondNumI;
			// Get the product
			productCluster = (NECluster *) network->getCompound(productName,
					compositionSizes);
		} else {
			// Just get the product if it is a single-species
			productCluster = (NECluster *) network->get(productName,
					productSize);
		}

		// React if the product exists in the network
		if (productCluster && (diffusionFactor > 0.0 || secondCluster->diffusionFactor > 0.0)) {
			// Setup the connectivity array for the second reactant
			setReactionConnectivity(secondCluster->id);
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondCluster, 0.0);
			// Push the product into the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);

			// Add it again if it is combining with itself
			if (secondCluster->id == id) combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void NECluster::replaceInCompound(std::vector<IReactant *> & reactants,
		const std::string& oldComponentName) {
	// Local Declarations
	std::map<std::string, int> secondReactantComp, productReactantComp;
	int numReactants = reactants.size();
	std::vector<int> productCompositionVector { 0, 0, 0 };
	NECluster *secondReactant, *productReactant;

	// Loop over all of the extra reactants in this reaction and handle the replacement
	for (int i = 0; i < numReactants; i++) {
		// Get the second reactant and its composition
		secondReactant = (NECluster *) reactants[i];
		secondReactantComp = secondReactant->getComposition();
		// Create the composition vector
		productReactantComp = secondReactantComp;
		// Updated the modified components
		productReactantComp[oldComponentName] =
				secondReactantComp[oldComponentName] - size;
		// Create the composition vector -- FIXME! This should be general!
		productCompositionVector = {productReactantComp[xeType],
				productReactantComp[vType],
				productReactantComp[iType] };
		// Get the product of the same type as the second reactant
		productReactant = (NECluster *) network->getCompound(secondReactant->typeName,
				productCompositionVector);
		// If the product exists, mark the proper reaction arrays and add it to the list
		if (productReactant) {
			// Setup the connectivity array for the second reactant
			setReactionConnectivity(secondReactant->id);
			// Creates the combining cluster
			// The reaction constant will be computed later and is set to 0.0 for now
			CombiningCluster combCluster(secondReactant, 0.0);
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);

			// Add it again if it is combining with itself
			if (secondReactant->id == id) combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void NECluster::fillVWithI(std::vector<IReactant *> & reactants) {
	// Local Declarations
	std::string productClusterName;
	int firstClusterSize = 0, secondClusterSize = 0, productClusterSize = 0,
			reactantVecSize = 0;
	NECluster *secondCluster, *productCluster;

	// Get the number of V or I in this cluster (the "first")
	firstClusterSize = size;
	// Look at all of the second clusters, either V or I, and determine
	// if a connection exists.
	reactantVecSize = reactants.size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second cluster its size
		secondCluster = (NECluster *) reactants[i];
		secondClusterSize = (secondCluster->size);
		// The only way this reaction is allowed is if the sizes are not equal.
		if (firstClusterSize != secondClusterSize) {
			// We have to switch on cluster type to make sure that the annihilation
			// is computed correctly.
			if (typeName == vType) {
				// Compute the product size and set the product name for the case
				// where I is the second cluster
				if (secondClusterSize > firstClusterSize) {
					productClusterSize = secondClusterSize - firstClusterSize;
					productClusterName = iType;
				} else if (secondClusterSize < firstClusterSize) {
					productClusterSize = firstClusterSize - secondClusterSize;
					productClusterName = vType;
				}
			} else if (typeName == iType) {
				// Compute the product size and set the product name for the case
				// where V is the second cluster
				if (firstClusterSize > secondClusterSize) {
					productClusterSize = firstClusterSize - secondClusterSize;
					productClusterName = iType;
				} else if (firstClusterSize < secondClusterSize) {
					productClusterSize = secondClusterSize - firstClusterSize;
					productClusterName = vType;
				}
			}
			// Get the product
			productCluster = (NECluster *) network->get(
					productClusterName, productClusterSize);
			// Only deal with this reaction if the product exists. Otherwise the
			// whole reaction is forbidden.
			if (productCluster) {
				// Setup the connectivity array to handle the second reactant
				setReactionConnectivity(secondCluster->id);
				// Creates the combining cluster
				// The reaction constant will be computed later and is set to 0.0 for now
				CombiningCluster combCluster(secondCluster, 0.0);
				// Push the second cluster onto the list of clusters that combine
				// with this one
				combiningReactants.push_back(combCluster);

				// Add it again if it is combining with itself
				if (secondCluster->id == id) combiningReactants.push_back(combCluster);
			}
		}
	}

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

std::vector<int> NECluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network->size());
}

std::vector<int> NECluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network->size());
}

void NECluster::resetConnectivities() {
	// Doesn't do anything because we are not using effective reactions anymore

	return;
}

void NECluster::setReactionNetwork(
		const std::shared_ptr<IReactionNetwork> reactionNetwork) {
	// Call the superclass's method to actually set the reference
	Reactant::setReactionNetwork(reactionNetwork);

	// Extract properties from the network
	auto properties = network->getProperties();

	// Get the enabled reaction type flags
	bool reactionsEnabled = (properties["reactionsEnabled"] == "true");
	bool dissociationsEnabled = (properties["dissociationsEnabled"] == "true");

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Get the index/id of this cluster in the reaction network.
	thisNetworkIndex = id - 1;

	// ----- Handle the connectivity for NEClusters -----

	// Generate the reactant and dissociation connectivity arrays.
	// This only must be done once since the arrays are stored as
	// member attributes. Only perform these tasks if the reaction
	// types are enabled.
	if (reactionsEnabled) {
		createReactionConnectivity();
	}
	if (dissociationsEnabled) {
		createDissociationConnectivity();
	}

	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	return;
}

double NECluster::getTotalFlux() {
	// Initialize the fluxes
	double prodFlux = 0.0, combFlux = 0.0, dissFlux = 0.0, emissFlux = 0.0;

	// Get the fluxes
	prodFlux = getProductionFlux();
	dissFlux = getDissociationFlux();

	// Don't compute the combination and emission flux if the
	// concentration is 0.0 because they are proportional to it
	if (!xolotlCore::equal(concentration, 0.0)) {
		combFlux = getCombinationFlux();
		emissFlux = getEmissionFlux();
	}

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double NECluster::getDissociationFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;
	NECluster *dissociatingCluster;

	// Set the total number of reactants that dissociate to form this one
	nPairs = dissociatingPairs.size();
	// Loop over all dissociating clusters that form this cluster
	for (int j = 0; j < nPairs; j++) {
		// Get the dissociating cluster
		dissociatingCluster = dissociatingPairs[j].first;
		// Calculate the Dissociation flux
		flux += dissociatingPairs[j].kConstant
				* dissociatingCluster->concentration;
	}

	// Return the flux
	return flux;
}

double NECluster::getEmissionFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;

	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop over all the pairs
	for (int i = 0; i < nPairs; i++) {
		// Update the flux
		flux += emissionPairs[i].kConstant;
	}

	return flux * concentration;
}

double NECluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	int nPairs = 0;
	NECluster *firstReactant, *secondReactant;

	// Set the total number of reacting pairs
	nPairs = reactingPairs.size();
	// Loop over all the reacting pairs
	for (int i = 0; i < nPairs; i++) {
		// Get the two reacting clusters
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Update the flux
		flux += reactingPairs[i].kConstant * firstReactant->concentration
				* secondReactant->concentration;
	}

	// Return the production flux
	return flux;
}

double NECluster::getCombinationFlux() const {
	// Local declarations
	double flux = 0.0;
	int nReactants = 0;
	NECluster *combiningCluster;

	// Set the total number of reactants that combine to form this one
	nReactants = combiningReactants.size();
	// Loop over all possible clusters
	for (int j = 0; j < nReactants; j++) {
		// Get the cluster that combines with this one
		combiningCluster = combiningReactants[j].combining;
		// Calculate Second term of production flux
		flux += combiningReactants[j].kConstant
				* combiningCluster->concentration;
	}

	return flux * concentration;
}

std::vector<double> NECluster::getPartialDerivatives() const {
	// Local Declarations
	std::vector<double> partials(network->size(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return partials;
}

void NECluster::getPartialDerivatives(std::vector<double> & partials) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void NECluster::getProductionPartialDerivatives(std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, index = 0;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	numReactants = reactingPairs.size();
	for (int i = 0; i < numReactants; i++) {
		// Compute the contribution from the first part of the reacting pair
		index = reactingPairs[i].first->id - 1;
		partials[index] += reactingPairs[i].kConstant
				* reactingPairs[i].second->concentration;
		// Compute the contribution from the second part of the reacting pair
		index = reactingPairs[i].second->id - 1;
		partials[index] += reactingPairs[i].kConstant
				* reactingPairs[i].first->concentration;
	}

	return;
}

void NECluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, otherIndex = 0;
	NECluster *cluster;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	numReactants = combiningReactants.size();
	for (int i = 0; i < numReactants; i++) {
		cluster = (NECluster *) combiningReactants[i].combining;
		// Get the index of cluster
		otherIndex = cluster->id - 1;
		// Remember that the flux due to combinations is OUTGOING (-=)!
		// Compute the contribution from this cluster
		partials[thisNetworkIndex] -= combiningReactants[i].kConstant
				* cluster->concentration;
		// Compute the contribution from the combining cluster
		partials[otherIndex] -= combiningReactants[i].kConstant
				* concentration;
	}

	return;
}

void NECluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0;
	NECluster *cluster;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	numPairs = dissociatingPairs.size();
	for (int i = 0; i < numPairs; i++) {
		// Get the dissociating cluster
		cluster = dissociatingPairs[i].first;
		index = cluster->id - 1;
		partials[index] += dissociatingPairs[i].kConstant;
	}

	return;
}

void NECluster::getEmissionPartialDerivatives(std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	numPairs = emissionPairs.size();
	for (int i = 0; i < numPairs; i++) {
		// Modify the partial derivative. Remember that the flux
		// due to emission is OUTGOING (-=)!
		index = id - 1;
		partials[index] -= emissionPairs[i].kConstant;
	}

	return;
}

void NECluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void NECluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double NECluster::getLeftSideRate() const {
	// Initialize the rate and the cluster pointer
	double totalRate = 0.0;
	NECluster *cluster;

	// Loop on the combining reactants
	for (int i = 0; i < combiningReactants.size(); i++) {
		cluster = (NECluster *) combiningReactants[i].combining;
		// Add the rate to the total rate
		totalRate += combiningReactants[i].kConstant
				* cluster->concentration;
	}

	// Loop on the emission pairs
	for (int i = 0; i < emissionPairs.size(); i++) {
		// Add the rate to the total rate
		totalRate += emissionPairs[i].kConstant;
	}

	return totalRate;
}

std::vector<int> NECluster::getConnectivity() const {
	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int)connectivityLength) {
		throw std::string("The reaction vector is an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int)connectivityLength) {
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

void NECluster::computeRateConstants() {
	// Local declarations
	NECluster *firstReactant, *secondReactant, *combiningReactant,
		*dissociatingCluster, *otherEmittedCluster, *firstCluster,
		*secondCluster;
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Compute the reaction constant associated to the reacting pairs
	// Set the total number of reacting pairs
	int nPairs = reactingPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		// Get the reactants
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*firstReactant,
				*secondReactant);
		// Set it in the pair
		reactingPairs[i].kConstant = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate) biggestProductionRate = rate;
	}

	// Compute the reaction constant associated to the combining reactants
	// Set the total number of combining reactants
	int nReactants = combiningReactants.size();
	// Loop on them
	for (int i = 0; i < nReactants; i++) {
		// Get the reactants
		combiningReactant = combiningReactants[i].combining;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*this, *combiningReactant);
		// Set it in the combining cluster
		combiningReactants[i].kConstant = rate;
	}

	// Compute the dissociation constant associated to the dissociating clusters
	// Set the total number of dissociating clusters
	nPairs = dissociatingPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		dissociatingCluster = dissociatingPairs[i].first;
		// The second element of the pair is the cluster that is also
		// emitted by the dissociation
		otherEmittedCluster = dissociatingPairs[i].second;
		// Compute the dissociation constant
		// The order of the cluster is important here because of the binding
		// energy used in the computation. It is taken from the type of the first cluster
		// which must be the single one
		if (size == 1) {
			// "this" is the single size one
			rate = calculateDissociationConstant(*dissociatingCluster, *this,
					*otherEmittedCluster);
		} else {
			// otherEmittedCluster is the single size one
			rate = calculateDissociationConstant(*dissociatingCluster,
					*otherEmittedCluster, *this);

		}
		// Set it in the pair
		dissociatingPairs[i].kConstant = rate;
	}

	// Compute the dissociation constant associated to the emission of pairs of clusters
	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop on them
	for (int i = 0; i < nPairs; i++) {
		firstCluster = emissionPairs[i].first;
		secondCluster = emissionPairs[i].second;
		// Compute the dissociation rate
		rate = calculateDissociationConstant(*this, *firstCluster,
				*secondCluster);
		// Set it in the pair
		emissionPairs[i].kConstant = rate;
	}

	// Shrink the arrays to save some space
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Set the biggest rate to the biggest production rate
	biggestRate = biggestProductionRate;

	return;
}
