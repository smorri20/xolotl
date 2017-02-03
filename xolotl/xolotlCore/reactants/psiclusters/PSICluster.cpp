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

void PSICluster::createReactionConnectivity() {
	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);

	// This cluster is always X_a

	// Initial declarations
	int firstSize = 0, secondSize = 0;

	// Single species clustering producing this cluster
	// X_(a-i) + X_i --> X_a
	for (firstSize = 1; firstSize <= (int) size / 2; firstSize++) {
		// Set the size of the second reactant
		secondSize = size - firstSize;
		// Get the first and second reactants for the reaction
		auto firstReactant = (PSICluster *) network->get(typeName, firstSize);
		auto secondReactant = (PSICluster *) network->get(typeName, secondSize);
		// Create a ReactingPair with the two reactants
		if (firstReactant && secondReactant
				&& (firstReactant->diffusionFactor > 0.0
						|| secondReactant->diffusionFactor > 0.0)) {
			// Create a production reaction
			auto reaction = std::make_shared<ProductionReaction>(firstReactant,
					secondReactant);
			// Add it to the network
			reaction = network->addProductionReaction(reaction);
			// Create a cluster pair
			ClusterPair pair(firstReactant, secondReactant, reaction.get());
			// Add the pair to the list
			reactingPairs.push_back(pair);
			// Setup the connectivity array
			setReactionConnectivity(firstReactant->id);
			setReactionConnectivity(secondReactant->id);
		}
	}

	// Single species clustering
	// X_a + X_b --> X_(a+b)
	auto reactants = network->getAll(typeName);
	// combineClusters handles everything for this type of reaction
	PSICluster::combineClusters(reactants, typeName);

	return;
}

void PSICluster::createDissociationConnectivity() {
	// This cluster is always X_a

	// Single species dissociation
	// X_a --> X_(a-1) + X
	auto smallerReactant = (PSICluster *) network->get(typeName, size - 1);
	auto singleCluster = (PSICluster *) network->get(typeName, 1);
	emitClusters(singleCluster, smallerReactant);

	// Single species dissociation producing this cluster
	// X_(a+1) --> X_a + X
	auto biggerReactant = (PSICluster *) network->get(typeName, size + 1);
	dissociateCluster(biggerReactant, singleCluster);

	// Specific case for the single size cluster
	// for a = 1
	if (size == 1) {
		// all the cluster of the same type dissociate into it
		auto allSameTypeReactants = network->getAll(typeName);
		for (unsigned int i = 0; i < allSameTypeReactants.size(); i++) {
			auto cluster = (PSICluster *) allSameTypeReactants[i];
			// X_1 cannot dissociate and X_2 --> X + X was already
			// counted in the previous step
			if (cluster->size < 3)
				continue;
			// X_b is the dissociating one, X_(b-a) is the one
			// that is also emitted during the dissociation
			smallerReactant = (PSICluster *) network->get(typeName,
					cluster->size - 1);
			dissociateCluster(cluster, smallerReactant);
		}
	}

	return;
}

void PSICluster::dissociateCluster(PSICluster * dissociatingCluster,
		PSICluster * emittedCluster) {
	// Test if the dissociatingCluster and the emittedCluster exist
	if (dissociatingCluster && emittedCluster
			&& (diffusionFactor > 0.0 || emittedCluster->diffusionFactor > 0.0)) {
		// Create a dissociation reaction
		auto reaction = std::make_shared<DissociationReaction>(
				dissociatingCluster, this, emittedCluster);
		// Add it to the network
		reaction = network->addDissociationReaction(reaction);
		// Create the pair of them where it is important that the
		// dissociating cluster is the first one
		ClusterPair pair(dissociatingCluster, emittedCluster, reaction.get());
		// Add the pair to the dissociating pair vector
		// The connectivity is handled in emitCluster
		dissociatingPairs.push_back(pair);

		// Take care of the connectivity
		setDissociationConnectivity(dissociatingCluster->id);

		// Add it to the list again if it the same as the other emitted cluster
		if (id == emittedCluster->id)
			dissociatingPairs.push_back(pair);
	}

	return;
}

void PSICluster::emitClusters(PSICluster * firstEmittedCluster,
		PSICluster * secondEmittedCluster) {
	// Test if the emitted clusters exist
	if (firstEmittedCluster && secondEmittedCluster
			&& (firstEmittedCluster->diffusionFactor > 0.0
					|| secondEmittedCluster->diffusionFactor > 0.0)) {
		// Connect this cluster to itself since any reaction will affect it
		setDissociationConnectivity(id);

		// Create a dissociation reaction
		auto reaction = std::make_shared<DissociationReaction>(this,
				firstEmittedCluster, secondEmittedCluster);
		// Add it to the network
		reaction = network->addDissociationReaction(reaction);

		// Add the pair of emitted clusters to the vector of emissionPairs
		// The first cluster is the size one one
		ClusterPair pair(firstEmittedCluster, secondEmittedCluster,
				reaction.get());
		emissionPairs.push_back(pair);
	}

	return;
}

void PSICluster::combineClusters(std::vector<IReactant *> & reactants,
		const std::string& productName) {
	// Initial declarations
	std::map<std::string, int> myComposition = getComposition(),
			secondComposition;
	int numHe = 0, numV = 0, numI = 0, secondNumHe = 0, secondNumV = 0,
			secondNumI = 0, productSize = 0;
	std::vector<int> compositionSizes { 0, 0, 0 };
	PSICluster *productCluster = nullptr, *secondCluster = nullptr;
	// Setup the composition variables for this cluster
	numHe = myComposition[heType];
	numV = myComposition[vType];
	numI = myComposition[iType];

	int reactantVecSize = reactants.size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second reactant, its composition and its index
		secondCluster = (PSICluster *) reactants[i];
		secondComposition = secondCluster->getComposition();
		secondNumHe = secondComposition[heType];
		secondNumV = secondComposition[vType];
		secondNumI = secondComposition[iType];
		// Compute the product size
		productSize = size + secondCluster->size;
		// Get and handle product for compounds
		if (productName == heVType || productName == heIType) {
			// Modify the composition vector
			compositionSizes[0] = numHe + secondNumHe;
			compositionSizes[1] = numV + secondNumV;
			compositionSizes[2] = numI + secondNumI;
			// Get the product
			productCluster = (PSICluster *) network->getCompound(productName,
					compositionSizes);
		} else {
			// Just get the product if it is a single-species
			productCluster = (PSICluster *) network->get(productName,
					productSize);
		}

		// React if the product exists in the network
		if (productCluster
				&& (diffusionFactor > 0.0
						|| secondCluster->diffusionFactor > 0.0)) {
			// Setup the connectivity array for the second reactant
			setReactionConnectivity(secondCluster->id);
			// Create the corresponding production reaction
			auto reaction = std::make_shared<ProductionReaction>(this,
					secondCluster);
			// Add it to the network
			reaction = network->addProductionReaction(reaction);
			// Creates the combining cluster
			CombiningCluster combCluster(secondCluster, reaction.get());
			// Push the product into the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);

			// Add it again if it is combining with itself
			if (secondCluster->id == id)
				combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void PSICluster::replaceInCompound(std::vector<IReactant *> & reactants,
		const std::string& oldComponentName) {
	// Local Declarations
	std::map<std::string, int> secondReactantComp, productReactantComp;
	int numReactants = reactants.size();
	std::vector<int> productCompositionVector { 0, 0, 0 };
	PSICluster *secondReactant = nullptr, *productReactant = nullptr;

	// Loop over all of the extra reactants in this reaction and handle the replacement
	for (int i = 0; i < numReactants; i++) {
		// Get the second reactant and its composition
		secondReactant = (PSICluster *) reactants[i];
		secondReactantComp = secondReactant->getComposition();
		// Create the composition vector
		productReactantComp = secondReactantComp;
		// Updated the modified components
		productReactantComp[oldComponentName] =
				secondReactantComp[oldComponentName] - size;
		// Create the composition vector -- FIXME! This should be general!
		productCompositionVector = {productReactantComp[heType],
			productReactantComp[vType],
			productReactantComp[iType]};
		// Get the product of the same type as the second reactant
		productReactant = (PSICluster *) network->getCompound(
				secondReactant->typeName, productCompositionVector);
		// If the product exists, mark the proper reaction arrays and add it to the list
		if (productReactant && (diffusionFactor > 0.0
				|| secondReactant->diffusionFactor > 0.0)) {
			// Setup the connectivity array for the second reactant
			setReactionConnectivity(secondReactant->id);
			// Create the corresponding production reaction
			auto reaction = std::make_shared<ProductionReaction>(this,
					secondReactant);
			// Add it to the network
			reaction = network->addProductionReaction(reaction);
			// Creates the combining cluster
			CombiningCluster combCluster(secondReactant, reaction.get());
			// Push the product onto the list of clusters that combine with this one
			combiningReactants.push_back(combCluster);

			// Add it again if it is combining with itself
			if (secondReactant->id == id)
				combiningReactants.push_back(combCluster);
		}
	}

	return;
}

void PSICluster::fillVWithI(std::vector<IReactant *> & reactants) {
	// Local Declarations
	std::string productClusterName;
	int firstClusterSize = 0, secondClusterSize = 0, productClusterSize = 0,
			reactantVecSize = 0;
	PSICluster *secondCluster = nullptr, *productCluster = nullptr;

	// Get the number of V or I in this cluster (the "first")
	firstClusterSize = size;
	// Look at all of the second clusters, either V or I, and determine
	// if a connection exists.
	reactantVecSize = reactants.size();
	for (int i = 0; i < reactantVecSize; i++) {
		// Get the second cluster its size
		secondCluster = (PSICluster *) reactants[i];
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
			productCluster = (PSICluster *) network->get(productClusterName,
					productClusterSize);
			// Only deal with this reaction if the product exists. Otherwise the
			// whole reaction is forbidden.
			if (productCluster && (diffusionFactor > 0.0
					|| secondCluster->diffusionFactor > 0.0)) {
				// Setup the connectivity array to handle the second reactant
				setReactionConnectivity(secondCluster->id);
				// Create the corresponding production reaction
				auto reaction = std::make_shared<ProductionReaction>(this,
						secondCluster);
				// Add it to the network
				reaction = network->addProductionReaction(reaction);
				// Creates the combining cluster
				CombiningCluster combCluster(secondCluster, reaction.get());
				// Push the product onto the list of clusters that combine with this one
				combiningReactants.push_back(combCluster);

				// Add it again if it is combining with itself
				if (secondCluster->id == id)
					combiningReactants.push_back(combCluster);
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
	for (auto it = reactingPairs.begin(); it != reactingPairs.end();
			++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->id);
		setReactionConnectivity((*it).second->id);
		setReactionConnectivity((*it).first->heMomId);
		setReactionConnectivity((*it).second->heMomId);
		setReactionConnectivity((*it).first->vMomId);
		setReactionConnectivity((*it).second->vMomId);
	}

	// Loop on the effective combining reactants
	for (auto it = combiningReactants.begin();
			it != combiningReactants.end(); ++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).combining->id);
		setReactionConnectivity((*it).combining->heMomId);
		setReactionConnectivity((*it).combining->vMomId);
	}

	// Loop on the effective dissociating pairs
	for (auto it = dissociatingPairs.begin();
			it != dissociatingPairs.end(); ++it) {
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

	// Get the enabled reaction type flags
	bool reactionsEnabled = reactionNetwork->getReactionsEnabled();
	bool dissociationsEnabled = reactionNetwork->getDissociationsEnabled();

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// ----- Handle the connectivity for PSIClusters -----

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

double PSICluster::getHeMomentum() const {
	return 0.0;
}

double PSICluster::getVMomentum() const {
	return 0.0;
}

double PSICluster::getTotalFlux() {
	// Get the fluxes
	double prodFlux = getProductionFlux();
	double dissFlux = getDissociationFlux();
	double combFlux = getCombinationFlux();
	double emissFlux = getEmissionFlux();

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double PSICluster::getDissociationFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;
	PSICluster *dissociatingCluster = nullptr;

	// Set the total number of reactants that dissociate to form this one
	nPairs = dissociatingPairs.size();
	// Loop over all dissociating clusters that form this cluster
	for (int j = 0; j < nPairs; j++) {
		// Get the dissociating cluster
		dissociatingCluster = dissociatingPairs[j].first;
		// Calculate the Dissociation flux
		flux += *(dissociatingPairs[j].kConstant)
				* dissociatingCluster->getConcentration(
						dissociatingPairs[j].firstHeDistance,
						dissociatingPairs[j].firstVDistance);
	}

	// Return the flux
	return flux;
}

double PSICluster::getEmissionFlux() const {
	// Initial declarations
	int nPairs = 0;
	double flux = 0.0;

	// Set the total number of emission pairs
	nPairs = emissionPairs.size();
	// Loop over all the pairs
	for (int i = 0; i < nPairs; i++) {
		// Update the flux
		flux += *(emissionPairs[i].kConstant);
	}

	return flux * concentration;
}

double PSICluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	int nPairs = 0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Set the total number of reacting pairs
	nPairs = reactingPairs.size();
	// Loop over all the reacting pairs
	for (int i = 0; i < nPairs; i++) {
		// Get the two reacting clusters
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Update the flux
		flux += *(reactingPairs[i].kConstant)
				* firstReactant->getConcentration(
						reactingPairs[i].firstHeDistance,
						reactingPairs[i].firstVDistance)
				* secondReactant->getConcentration(
						reactingPairs[i].secondHeDistance,
						reactingPairs[i].secondVDistance);
	}

	// Return the production flux
	return flux;
}

double PSICluster::getCombinationFlux() const {
	// Local declarations
	double flux = 0.0;
	int nReactants = 0;
	PSICluster *combiningCluster = nullptr;

	// Set the total number of reactants that combine to form this one
	nReactants = combiningReactants.size();
	// Loop over all possible clusters
	for (int j = 0; j < nReactants; j++) {
		// Get the cluster that combines with this one
		combiningCluster = combiningReactants[j].combining;
		// Calculate the combination flux
		flux += *(combiningReactants[j].kConstant)
				* combiningCluster->getConcentration(
						combiningReactants[j].heDistance,
						combiningReactants[j].vDistance);
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
	int numReactants = 0, index = 0;
	double value = 0.0;

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
		value = *(reactingPairs[i].kConstant)
				* reactingPairs[i].second->getConcentration(
						reactingPairs[i].secondHeDistance,
						reactingPairs[i].secondVDistance);
		index = reactingPairs[i].first->id - 1;
		partials[index] += value;
		index = reactingPairs[i].first->heMomId - 1;
		partials[index] += value * reactingPairs[i].firstHeDistance;
		index = reactingPairs[i].first->vMomId - 1;
		partials[index] += value * reactingPairs[i].firstVDistance;
		// Compute the contribution from the second part of the reacting pair
		value = *(reactingPairs[i].kConstant)
				* reactingPairs[i].first->getConcentration(
						reactingPairs[i].firstHeDistance,
						reactingPairs[i].firstVDistance);
		index = reactingPairs[i].second->id - 1;
		partials[index] += value;
		index = reactingPairs[i].second->heMomId - 1;
		partials[index] += value * reactingPairs[i].secondHeDistance;
		index = reactingPairs[i].second->vMomId - 1;
		partials[index] += value * reactingPairs[i].secondVDistance;
	}

	return;
}

void PSICluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, otherIndex = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	numReactants = combiningReactants.size();
	for (int i = 0; i < numReactants; i++) {
		cluster = (PSICluster *) combiningReactants[i].combining;
		// Remember that the flux due to combinations is OUTGOING (-=)!
		// Compute the contribution from this cluster
		partials[id - 1] -= *(combiningReactants[i].kConstant)
				* cluster->getConcentration(
						combiningReactants[i].heDistance,
						combiningReactants[i].vDistance);
		// Compute the contribution from the combining cluster
		value = *(combiningReactants[i].kConstant) * concentration;
		otherIndex = cluster->id - 1;
		partials[otherIndex] -= value;
		otherIndex = cluster->heMomId - 1;
		partials[otherIndex] -= value * combiningReactants[i].heDistance;
		otherIndex = cluster->vMomId - 1;
		partials[otherIndex] -= value * combiningReactants[i].vDistance;
	}

	return;
}

void PSICluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

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
		value = *(dissociatingPairs[i].kConstant);
		index = cluster->id - 1;
		partials[index] += value;
		index = cluster->heMomId - 1;
		partials[index] += value * dissociatingPairs[i].firstHeDistance;
		index = cluster->vMomId - 1;
		partials[index] += value * dissociatingPairs[i].firstVDistance;
	}

	return;
}

void PSICluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
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
		partials[index] -= *(emissionPairs[i].kConstant);
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
		totalRate += *(combiningReactants[i].kConstant)
				* cluster->concentration;
	}

	// Loop on the emission pairs
	for (int i = 0; i < emissionPairs.size(); i++) {
		// Add the rate to the total rate
		totalRate += *(emissionPairs[i].kConstant);
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
