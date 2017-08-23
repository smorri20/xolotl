// Includes
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <MathUtils.h>

using namespace xolotlCore;

/**
 * The helium momentum partials.
 */
std::vector<double> heMomentumPartials;

/**
 * The vacancy momentum partials.
 */
std::vector<double> vMomentumPartials;

PSISuperCluster::PSISuperCluster(double numHe, double numV, int nTot,
		int heWidth, int vWidth,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry), numHe(numHe), numV(numV), nTot(nTot), lowerHe(0), upperHe(
				0), lowerV(0), upperV(0), l0(0.0), l1He(0.0), l1V(0.0), dispersionHe(
				0.0), dispersionV(0.0), heMomentumFlux(0.0), vMomentumFlux(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	compositionMap[Species::He] = (int) numHe;
	compositionMap[Species::V] = (int) numV;

	// Set the width
	sectionHeWidth = heWidth;
	sectionVWidth = vWidth;

	// Set the formation energy
	formationEnergy = 0.0; // It is set to 0.0 because we do not want the super clusters to undergo dissociation

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	type = Species::PSISuper;

	return;
}

PSISuperCluster::PSISuperCluster(PSISuperCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;
	nTot = other.nTot;
	sectionHeWidth = other.sectionHeWidth;
	sectionVWidth = other.sectionVWidth;
	l0 = other.l0;
	l1He = other.l1He;
	l1V = other.l1V;
	lowerHe = other.lowerHe;
	upperHe = other.upperHe;
	lowerV = other.lowerV;
	upperV = other.upperV;
	dispersionHe = other.dispersionHe;
	dispersionV = other.dispersionV;
	effReactingList = other.effReactingList;
	effCombiningList = other.effCombiningList;
	effDissociatingList = other.effDissociatingList;
	effEmissionList = other.effEmissionList;
	heMomentumFlux = other.heMomentumFlux;
	vMomentumFlux = other.vMomentumFlux;

	return;
}

void PSISuperCluster::createProduction(
		std::shared_ptr<ProductionReaction> reaction, int a, int b, int c,
		int d) {
	// Check if the reaction was already added
	std::forward_list<SuperClusterProductionPair>::iterator it;
	for (it = effReactingList.begin(); it != effReactingList.end(); it++) {
		if (reaction->first == (*it).first
				&& reaction->second == (*it).second) {
			break;
		}
	}
	if (it == effReactingList.end()) {
		// It was not already in so add it
		// Add the production reaction to the network
		reaction = network.addProductionReaction(reaction);
		// Create a new SuperClusterProductionPair
		SuperClusterProductionPair superPair((PSICluster *) reaction->first,
				(PSICluster *) reaction->second, reaction.get());
		// Add it
		effReactingList.push_front(superPair);
		it = effReactingList.begin();
	}

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
			secondVDistance = 0.0;
	if (reaction->first->getType() == Species::PSISuper) {
		auto super = (PSICluster *) reaction->first;
		firstHeDistance = super->getHeDistance(c);
		firstVDistance = super->getVDistance(d);
	}
	if (reaction->second->getType() == Species::PSISuper) {
		auto super = (PSICluster *) reaction->second;
		secondHeDistance = super->getHeDistance(c);
		secondVDistance = super->getVDistance(d);
	}
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vFactor = (double) (b - numV) / dispersionV;
	// First is A, second is B, in A + B -> this
	(*it).a000 += 1.0;
	(*it).a001 += heFactor;
	(*it).a002 += vFactor;
	(*it).a100 += firstHeDistance;
	(*it).a101 += firstHeDistance * heFactor;
	(*it).a102 += firstHeDistance * vFactor;
	(*it).a200 += firstVDistance;
	(*it).a201 += firstVDistance * heFactor;
	(*it).a202 += firstVDistance * vFactor;
	(*it).a010 += secondHeDistance;
	(*it).a011 += secondHeDistance * heFactor;
	(*it).a012 += secondHeDistance * vFactor;
	(*it).a020 += secondVDistance;
	(*it).a021 += secondVDistance * heFactor;
	(*it).a022 += secondVDistance * vFactor;
	(*it).a110 += firstHeDistance * secondHeDistance;
	(*it).a111 += firstHeDistance * secondHeDistance * heFactor;
	(*it).a112 += firstHeDistance * secondHeDistance * vFactor;
	(*it).a120 += firstHeDistance * secondVDistance;
	(*it).a121 += firstHeDistance * secondVDistance * heFactor;
	(*it).a122 += firstHeDistance * secondVDistance * vFactor;
	(*it).a210 += firstVDistance * secondHeDistance;
	(*it).a211 += firstVDistance * secondHeDistance * heFactor;
	(*it).a212 += firstVDistance * secondHeDistance * vFactor;
	(*it).a220 += firstVDistance * secondVDistance;
	(*it).a221 += firstVDistance * secondVDistance * heFactor;
	(*it).a222 += firstVDistance * secondVDistance * vFactor;

	return;
}

void PSISuperCluster::createCombination(
		std::shared_ptr<ProductionReaction> reaction, int a, int b) {
	setReactionConnectivity(id);
	// Look for the other cluster
	IReactant * secondCluster;
	if (reaction->first->getId() == id)
		secondCluster = reaction->second;
	else
		secondCluster = reaction->first;

	// Check if the reaction was already added
	std::forward_list<SuperClusterProductionPair>::iterator it;
	for (it = effCombiningList.begin(); it != effCombiningList.end(); it++) {
		if (secondCluster == (*it).first) {
			break;
		}
	}
	if (it == effCombiningList.end()) {
		// It was not already in so add it
		// Create the corresponding production reaction
		auto newReaction = std::make_shared<ProductionReaction>(this,
				secondCluster);
		// Add it to the network
		newReaction = network.addProductionReaction(newReaction);
		// Create a new SuperClusterProductionPair
		SuperClusterProductionPair superPair((PSICluster *) secondCluster,
				nullptr, newReaction.get());
		// Add it
		effCombiningList.push_front(superPair);
		it = effCombiningList.begin();
	}

	// Update the coefficients
	double heDistance = getHeDistance(a);
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vDistance = getVDistance(b);
	double vFactor = (double) (b - numV) / dispersionV;
	// This is A, itBis is B, in A + B -> C
	(*it).a000 += 1.0;
	(*it).a001 += heFactor;
	(*it).a002 += vFactor;
	(*it).a100 += heDistance;
	(*it).a101 += heDistance * heFactor;
	(*it).a102 += heDistance * vFactor;
	(*it).a200 += vDistance;
	(*it).a201 += vDistance * heFactor;
	(*it).a202 += vDistance * vFactor;

	return;
}

void PSISuperCluster::createDissociation(
		std::shared_ptr<DissociationReaction> reaction, int a, int b, int c,
		int d) {
	// Look for the other cluster
	IReactant * emittedCluster;
	if (reaction->first->getId() == id)
		emittedCluster = reaction->second;
	else
		emittedCluster = reaction->first;

	// Check if the reaction was already added
	std::forward_list<SuperClusterDissociationPair>::iterator it;
	for (it = effDissociatingList.begin(); it != effDissociatingList.end();
			it++) {
		if (reaction->dissociating == (*it).first
				&& emittedCluster == (*it).second) {
			break;
		}
	}
	if (it == effDissociatingList.end()) {
		// It was not already in so add it
		// Create a dissociation reaction
		auto newReaction = std::make_shared<DissociationReaction>(
				reaction->dissociating, this, emittedCluster);
		// Add it to the network
		newReaction = network.addDissociationReaction(newReaction);
		// Create a new SuperClusterDissociationPair
		SuperClusterDissociationPair superPair(
				(PSICluster *) reaction->dissociating,
				(PSICluster *) emittedCluster, newReaction.get());
		// Add it
		effDissociatingList.push_front(superPair);
		it = effDissociatingList.begin();
	}

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction->dissociating->getType() == Species::PSISuper) {
		auto super = (PSICluster *) reaction->dissociating;
		firstHeDistance = super->getHeDistance(a);
		firstVDistance = super->getVDistance(b);
	}
	double heFactor = (double) (c - numHe) / dispersionHe;
	double vFactor = (double) (d - numV) / dispersionV;

	// A is the dissociating cluster
	(*it).a00 += 1.0;
	(*it).a01 += heFactor;
	(*it).a02 += vFactor;
	(*it).a10 += firstHeDistance;
	(*it).a11 += firstHeDistance * heFactor;
	(*it).a12 += firstHeDistance * vFactor;
	(*it).a20 += firstVDistance;
	(*it).a21 += firstVDistance * heFactor;
	(*it).a22 += firstVDistance * vFactor;

	return;
}

void PSISuperCluster::createEmission(
		std::shared_ptr<DissociationReaction> reaction, int a, int b, int c,
		int d) {
	// Check if the reaction was already added
	std::forward_list<SuperClusterDissociationPair>::iterator it;
	for (it = effEmissionList.begin(); it != effEmissionList.end(); it++) {
		if (reaction->first == (*it).first
				&& reaction->second == (*it).second) {
			break;
		}
	}
	if (it == effEmissionList.end()) {
		// It was not already in so add it
		// Add the reaction to the network
		reaction = network.addDissociationReaction(reaction);
		// Create a new SuperClusterDissociationPair
		SuperClusterDissociationPair superPair((PSICluster *) reaction->first,
				(PSICluster *) reaction->second, reaction.get());
		// Add it
		effEmissionList.push_front(superPair);
		it = effEmissionList.begin();
	}

	// Update the coeeficients
	double heDistance = getHeDistance(a);
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vDistance = getVDistance(b);
	double vFactor = (double) (b - numV) / dispersionV;
	// A is the dissociating cluster
	(*it).a00 += 1.0;
	(*it).a01 += heFactor;
	(*it).a02 += vFactor;
	(*it).a10 += heDistance;
	(*it).a11 += heDistance * heFactor;
	(*it).a12 += heDistance * vFactor;
	(*it).a20 += vDistance;
	(*it).a21 += vDistance * heFactor;
	(*it).a22 += vDistance * vFactor;

	return;
}

void PSISuperCluster::updateFromNetwork() {

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

void PSISuperCluster::setHeVVector(std::vector<std::pair<int, int> > vec) {
	// Initialize the dispersion sum
	double nHeSquare = 0.0, nVSquare = 0.0;
	// Update the network map, compute the radius and dispersions
	for (auto it = vec.begin(); it != vec.end(); it++) {
		reactionRadius += xolotlCore::tungstenLatticeConstant
				* pow((3.0 * (double) ((*it).second)) / xolotlCore::pi,
						(1.0 / 3.0)) * 0.5 / (double) nTot;

		// Compute nSquare for the dispersion
		nHeSquare += (double) (*it).first * (*it).first;
		nVSquare += (double) (*it).second * (*it).second;
	}

	// Compute the dispersions
	if (sectionHeWidth == 1)
		dispersionHe = 1.0;
	else
		dispersionHe = 2.0 * (nHeSquare - (numHe * (double) nTot * numHe))
				/ ((double) (nTot * (sectionHeWidth - 1)));

	if (sectionVWidth == 1)
		dispersionV = 1.0;
	else
		dispersionV = 2.0 * (nVSquare - (numV * (double) nTot * numV))
				/ ((double) (nTot * (sectionVWidth - 1)));

	// Set the boundaries
	lowerHe = (int) (numHe - (double) sectionHeWidth / 2.0) + 1;
	upperHe = (int) (numHe - (double) sectionHeWidth / 2.0) + sectionHeWidth;
	lowerV = (int) (numV - (double) sectionVWidth / 2.0) + 1;
	upperV = (int) (numV - (double) sectionVWidth / 2.0) + sectionVWidth;

	return;
}

double PSISuperCluster::getTotalConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (int i = lowerHe; i <= upperHe; i++) {
		for (int j = lowerV; j <= upperV; j++) {
			// Compute the distances
			heDistance = getHeDistance(i);
			vDistance = getVDistance(j);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance);
		}
	}

	return conc;
}

double PSISuperCluster::getTotalHeliumConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (int i = lowerHe; i <= upperHe; i++) {
		for (int j = lowerV; j <= upperV; j++) {
			// Compute the distances
			heDistance = getHeDistance(i);
			vDistance = getVDistance(j);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) i;
		}
	}

	return conc;
}

double PSISuperCluster::getTotalVacancyConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (int i = lowerHe; i <= upperHe; i++) {
		for (int j = lowerV; j <= upperV; j++) {
			// Compute the distances
			heDistance = getHeDistance(i);
			vDistance = getVDistance(j);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) j;
		}
	}

	return conc;
}

void PSISuperCluster::resetConnectivities() {
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

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// The cluster is connecting to both clusters in the pair
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getHeMomentumId());
		setReactionConnectivity((*it).first->getVMomentumId());
		setReactionConnectivity((*it).second->getId());
		setReactionConnectivity((*it).second->getHeMomentumId());
		setReactionConnectivity((*it).second->getVMomentumId());
	}

	// Loop over all the combining pairs
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getHeMomentumId());
		setReactionConnectivity((*it).first->getVMomentumId());
	}

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// The cluster is connecting to the combining cluster
		setDissociationConnectivity((*it).first->getId());
		setDissociationConnectivity((*it).first->getHeMomentumId());
		setDissociationConnectivity((*it).first->getVMomentumId());
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the momentum
	int dof = network.getDOF();
	heMomentumPartials.resize(dof, 0.0);
	vMomentumPartials.resize(dof, 0.0);

	return;
}

double PSISuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;
	PSICluster *dissociatingCluster = nullptr;

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		dissociatingCluster = (*it).first;
		double l0A = dissociatingCluster->getConcentration(0.0, 0.0);
		double lHeA = dissociatingCluster->getHeMomentum();
		double lVA = dissociatingCluster->getVMomentum();
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value * ((*it).a00 * l0A + (*it).a10 * lHeA + (*it).a20 * lVA);
		// Compute the momentum fluxes
		heMomentumFlux += value
				* ((*it).a01 * l0A + (*it).a11 * lHeA + (*it).a21 * lVA);
		vMomentumFlux += value
				* ((*it).a02 * l0A + (*it).a12 * lHeA + (*it).a22 * lVA);
	}

	// Return the flux
	return flux;
}

double PSISuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value * ((*it).a00 * l0 + (*it).a10 * l1He + (*it).a20 * l1V);
		// Compute the momentum fluxes
		heMomentumFlux -= value
				* ((*it).a01 * l0 + (*it).a11 * l1He + (*it).a21 * l1V);
		vMomentumFlux -= value
				* ((*it).a02 * l0 + (*it).a12 * l1He + (*it).a22 * l1V);
	}

	return flux;
}

double PSISuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
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
		value = *((*it).kConstant) / (double) nTot;
		flux += value
				* ((*it).a000 * l0A * l0B + (*it).a010 * l0A * lHeB
						+ (*it).a020 * l0A * lVB + (*it).a100 * lHeA * l0B
						+ (*it).a110 * lHeA * lHeB + (*it).a120 * lHeA * lVB
						+ (*it).a200 * lVA * l0B + (*it).a210 * lVA * lHeB
						+ (*it).a220 * lVA * lVB);
		// Compute the momentum fluxes
		heMomentumFlux += value
				* ((*it).a001 * l0A * l0B + (*it).a011 * l0A * lHeB
						+ (*it).a021 * l0A * lVB + (*it).a101 * lHeA * l0B
						+ (*it).a111 * lHeA * lHeB + (*it).a121 * lHeA * lVB
						+ (*it).a201 * lVA * l0B + (*it).a211 * lVA * lHeB
						+ (*it).a221 * lVA * lVB);
		vMomentumFlux += value
				* ((*it).a002 * l0A * l0B + (*it).a012 * l0A * lHeB
						+ (*it).a022 * l0A * lVB + (*it).a102 * lHeA * l0B
						+ (*it).a112 * lHeA * lHeB + (*it).a122 * lHeA * lVB
						+ (*it).a202 * lVA * l0B + (*it).a212 * lVA * lHeB
						+ (*it).a222 * lVA * lVB);
	}

	// Return the production flux
	return flux;
}

double PSISuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	PSICluster *combiningCluster = nullptr;

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the combining cluster
		combiningCluster = (*it).first;
		double l0B = combiningCluster->getConcentration(0.0, 0.0);
		double lHeB = combiningCluster->getHeMomentum();
		double lVB = combiningCluster->getVMomentum();
		// Update the flux
		value = *((*it).kConstant) / (double) nTot;
		flux += value
				* ((*it).a000 * l0B * l0 + (*it).a100 * l0B * l1He
						+ (*it).a200 * l0B * l1V + (*it).a010 * lHeB * l0
						+ (*it).a110 * lHeB * l1He + (*it).a210 * lHeB * l1V
						+ (*it).a020 * lVB * l0 + (*it).a120 * lVB * l1He
						+ (*it).a220 * lVB * l1V);
		// Compute the momentum fluxes
		heMomentumFlux -= value
				* ((*it).a001 * l0B * l0 + (*it).a101 * l0B * l1He
						+ (*it).a201 * l0B * l1V + (*it).a011 * lHeB * l0
						+ (*it).a111 * lHeB * l1He + (*it).a211 * lHeB * l1V
						+ (*it).a021 * lVB * l0 + (*it).a121 * lVB * l1He
						+ (*it).a221 * lVB * l1V);
		vMomentumFlux -= value
				* ((*it).a002 * l0B * l0 + (*it).a102 * l0B * l1He
						+ (*it).a202 * l0B * l1V + (*it).a012 * lHeB * l0
						+ (*it).a112 * lHeB * l1He + (*it).a212 * lHeB * l1V
						+ (*it).a022 * lVB * l0 + (*it).a122 * lVB * l1He
						+ (*it).a222 * lVB * l1V);
	}

	return flux;
}

void PSISuperCluster::getPartialDerivatives(
		std::vector<double> & partials) const {
	// Reinitialize the momentum partial derivatives vector
	std::fill(heMomentumPartials.begin(), heMomentumPartials.end(), 0.0);
	std::fill(vMomentumPartials.begin(), vMomentumPartials.end(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSISuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	double value = 0.0;
	int index = 0;
	PSICluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
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
		value = *((*it).kConstant) / (double) nTot;
		index = firstReactant->getId() - 1;
		partials[index] += value
				* ((*it).a000 * l0B + (*it).a010 * lHeB + (*it).a020 * lVB);
		heMomentumPartials[index] += value
				* ((*it).a001 * l0B + (*it).a011 * lHeB + (*it).a021 * lVB);
		vMomentumPartials[index] += value
				* ((*it).a002 * l0B + (*it).a012 * lHeB + (*it).a022 * lVB);
		index = firstReactant->getHeMomentumId() - 1;
		partials[index] += value
				* ((*it).a100 * l0B + (*it).a110 * lHeB + (*it).a120 * lVB);
		heMomentumPartials[index] += value
				* ((*it).a101 * l0B + (*it).a111 * lHeB + (*it).a121 * lVB);
		vMomentumPartials[index] += value
				* ((*it).a102 * l0B + (*it).a112 * lHeB + (*it).a122 * lVB);
		index = firstReactant->getVMomentumId() - 1;
		partials[index] += value
				* ((*it).a200 * l0B + (*it).a210 * lHeB + (*it).a220 * lVB);
		heMomentumPartials[index] += value
				* ((*it).a201 * l0B + (*it).a211 * lHeB + (*it).a221 * lVB);
		vMomentumPartials[index] += value
				* ((*it).a202 * l0B + (*it).a212 * lHeB + (*it).a222 * lVB);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant->getId() - 1;
		partials[index] += value
				* ((*it).a000 * l0A + (*it).a100 * lHeA + (*it).a200 * lVA);
		heMomentumPartials[index] += value
				* ((*it).a001 * l0A + (*it).a101 * lHeA + (*it).a201 * lVA);
		vMomentumPartials[index] += value
				* ((*it).a002 * l0A + (*it).a102 * lHeA + (*it).a202 * lVA);
		index = secondReactant->getHeMomentumId() - 1;
		partials[index] += value
				* ((*it).a010 * l0A + (*it).a110 * lHeA + (*it).a210 * lVA);
		heMomentumPartials[index] += value
				* ((*it).a011 * l0A + (*it).a111 * lHeA + (*it).a211 * lVA);
		vMomentumPartials[index] += value
				* ((*it).a012 * l0A + (*it).a112 * lHeA + (*it).a212 * lVA);
		index = secondReactant->getVMomentumId() - 1;
		partials[index] += value
				* ((*it).a020 * l0A + (*it).a120 * lHeA + (*it).a220 * lVA);
		heMomentumPartials[index] += value
				* ((*it).a021 * l0A + (*it).a121 * lHeA + (*it).a221 * lVA);
		vMomentumPartials[index] += value
				* ((*it).a022 * l0A + (*it).a122 * lHeA + (*it).a222 * lVA);
	}

	return;
}

void PSISuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	PSICluster *cluster = nullptr;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end();
			++it) {
		// Get the combining clusters
		cluster = (*it).first;
		double l0B = cluster->getConcentration(0.0, 0.0);
		double lHeB = cluster->getHeMomentum();
		double lVB = cluster->getVMomentum();

		// Compute the contribution from the combining cluster
		value = *((*it).kConstant) / (double) nTot;
		index = cluster->getId() - 1;
		partials[index] -= value
				* ((*it).a000 * l0 + (*it).a100 * l1He + (*it).a200 * l1V);
		heMomentumPartials[index] -= value
				* ((*it).a001 * l0 + (*it).a101 * l1He + (*it).a201 * l1V);
		vMomentumPartials[index] -= value
				* ((*it).a002 * l0 + (*it).a102 * l1He + (*it).a202 * l1V);
		index = cluster->getHeMomentumId() - 1;
		partials[index] -= value
				* ((*it).a010 * l0 + (*it).a110 * l1He + (*it).a210 * l1V);
		heMomentumPartials[index] -= value
				* ((*it).a011 * l0 + (*it).a111 * l1He + (*it).a211 * l1V);
		vMomentumPartials[index] -= value
				* ((*it).a012 * l0 + (*it).a112 * l1He + (*it).a212 * l1V);
		index = cluster->getVMomentumId() - 1;
		partials[index] -= value
				* ((*it).a020 * l0 + (*it).a120 * l1He + (*it).a220 * l1V);
		heMomentumPartials[index] -= value
				* ((*it).a021 * l0 + (*it).a121 * l1He + (*it).a221 * l1V);
		vMomentumPartials[index] -= value
				* ((*it).a022 * l0 + (*it).a122 * l1He + (*it).a222 * l1V);
		// Compute the contribution from this cluster
		index = id - 1;
		partials[index] -= value
				* ((*it).a000 * l0B + (*it).a010 * lHeB + (*it).a020 * lVB);
		heMomentumPartials[index] -= value
				* ((*it).a001 * l0B + (*it).a011 * lHeB + (*it).a021 * lVB);
		vMomentumPartials[index] -= value
				* ((*it).a002 * l0B + (*it).a012 * lHeB + (*it).a022 * lVB);
		index = heMomId - 1;
		partials[index] -= value
				* ((*it).a100 * l0B + (*it).a110 * lHeB + (*it).a120 * lVB);
		heMomentumPartials[index] -= value
				* ((*it).a101 * l0B + (*it).a111 * lHeB + (*it).a121 * lVB);
		vMomentumPartials[index] -= value
				* ((*it).a102 * l0B + (*it).a112 * lHeB + (*it).a122 * lVB);
		index = vMomId - 1;
		partials[index] -= value
				* ((*it).a200 * l0B + (*it).a210 * lHeB + (*it).a220 * lVB);
		heMomentumPartials[index] -= value
				* ((*it).a201 * l0B + (*it).a211 * lHeB + (*it).a221 * lVB);
		vMomentumPartials[index] -= value
				* ((*it).a202 * l0B + (*it).a212 * lHeB + (*it).a222 * lVB);
	}

	return;
}

void PSISuperCluster::getDissociationPartialDerivatives(
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

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end();
			++it) {
		// Get the dissociating clusters
		cluster = (*it).first;
		// Compute the contribution from the dissociating cluster
		value = *((*it).kConstant) / (double) nTot;
		index = cluster->getId() - 1;
		partials[index] += value * ((*it).a00);
		heMomentumPartials[index] += value * ((*it).a01);
		vMomentumPartials[index] += value * ((*it).a02);
		index = cluster->getHeMomentumId() - 1;
		partials[index] += value * ((*it).a10);
		heMomentumPartials[index] += value * ((*it).a11);
		vMomentumPartials[index] += value * ((*it).a12);
		index = cluster->getVMomentumId() - 1;
		partials[index] += value * ((*it).a20);
		heMomentumPartials[index] += value * ((*it).a21);
		vMomentumPartials[index] += value * ((*it).a22);
	}

	return;
}

void PSISuperCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	double value = 0.0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Compute the contribution from the dissociating cluster
		value = *((*it).kConstant) / (double) nTot;
		index = id - 1;
		partials[index] -= value * ((*it).a00);
		heMomentumPartials[index] -= value * ((*it).a01);
		vMomentumPartials[index] -= value * ((*it).a02);
		index = heMomId - 1;
		partials[index] -= value * ((*it).a10);
		heMomentumPartials[index] -= value * ((*it).a11);
		vMomentumPartials[index] -= value * ((*it).a12);
		index = vMomId - 1;
		partials[index] -= value * ((*it).a20);
		heMomentumPartials[index] -= value * ((*it).a21);
		vMomentumPartials[index] -= value * ((*it).a22);
	}

	return;
}

void PSISuperCluster::getHeMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = heMomentumPartials[i];
	}

	return;
}

void PSISuperCluster::getVMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = vMomentumPartials[i];
	}

	return;
}
