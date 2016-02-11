// Includes
#include "SuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

SuperCluster::SuperCluster(double numHe, double numV, int nTot, int heWidth,
		int vWidth, double radius, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHe), numV(numV), nTot(nTot), l0(0.0), l1He(
				0.0), l1V(0.0), dispersionHe(0.0), dispersionV(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	compositionMap[heType] = (int) (numHe * (double) nTot);
	compositionMap[vType] = (int) (numV * (double) nTot);

	// Set the width
	sectionHeWidth = heWidth;
	sectionVWidth = vWidth;

	// Set the reaction radius and formation energy
	reactionRadius = radius;
	formationEnergy = 0.0;

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = "Super";

	return;
}

SuperCluster::SuperCluster(const SuperCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;
	nTot = other.nTot;
	sectionHeWidth = other.sectionHeWidth;
	sectionVWidth = other.sectionVWidth;
	l0 = other.l0;
	l1He = other.l1He;
	l1V = other.l1V;
	dispersionHe = other.dispersionHe;
	dispersionV = other.dispersionV;
	reactingMap = other.reactingMap;
	combiningMap = other.combiningMap;
	dissociatingMap = other.dissociatingMap;
	emissionMap = other.emissionMap;
	effReactingMap = other.effReactingMap;
	effCombiningMap = other.effCombiningMap;
	effDissociatingMap = other.effDissociatingMap;
	effEmissionMap = other.effEmissionMap;

	return;
}

std::shared_ptr<Reactant> SuperCluster::clone() {
	std::shared_ptr<Reactant> reactant(new SuperCluster(*this));

	return reactant;
}

double SuperCluster::getConcentration(double distHe, double distV) const {
	return l0 + (distHe * l1He) + (distV * l1V);
}

double SuperCluster::getTotalConcentration() const {
	// Initial declarations
	int heIndex = 0, vIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Check if this cluster exists
			if (effReactingMap.find(std::make_pair(heIndex, vIndex)) == effReactingMap.end())
				continue;

			// Compute the distances
			heDistance = (double) heIndex - numHe;
			vDistance = (double) vIndex - numV;

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance);
		}
	}

	return conc;
}

double SuperCluster::getTotalHeliumConcentration() const {
	// Initial declarations
	int heIndex = 0, vIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Check if this cluster exists
			if (effReactingMap.find(std::make_pair(heIndex, vIndex)) == effReactingMap.end())
				continue;

			// Compute the distances
			heDistance = (double) heIndex - numHe;
			vDistance = (double) vIndex - numV;

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) heIndex;
		}
	}

	return conc;
}

double SuperCluster::getHeDistance(int he) const {
	return (double) he - numHe;
}

double SuperCluster::getVDistance(int v) const {
	return (double) v - numV;
}

void SuperCluster::createReactionConnectivity() {
	// Aggregate the reacting pairs and combining reactants from the heVVector
	// Loop on the heVVector
	for (int i = 0; i < heVVector.size(); i++) {
		// Get the cluster composition
		auto comp = heVVector[i]->getComposition();
		// Get both production vectors
		auto react = heVVector[i]->reactingPairs;
		auto combi = heVVector[i]->combiningReactants;

		// Create the key to the map
		auto key = std::make_pair(comp[heType], comp[vType]);

		// Set them in the super cluster map
		reactingMap[key] = react;
		combiningMap[key] = combi;
	}

	return;
}

void SuperCluster::createDissociationConnectivity() {
	// Aggregate the dissociating and emission pairs from the heVVector
	// Loop on the heVVector
	for (int i = 0; i < heVVector.size(); i++) {
		// Get the cluster composition
		auto comp = heVVector[i]->getComposition();
		// Get both dissociation vectors
		auto disso = heVVector[i]->dissociatingPairs;
		auto emi = heVVector[i]->emissionPairs;

		// Create the key to the map
		auto key = std::make_pair(comp[heType], comp[vType]);

		// Set them in the super cluster map
		dissociatingMap[key] = disso;
		emissionMap[key] = emi;
	}

	return;
}

void SuperCluster::computeRateConstants() {
	// Local declarations
	PSICluster *firstReactant, *secondReactant, *combiningReactant,
			*dissociatingCluster, *otherEmittedCluster, *firstCluster,
			*secondCluster;
	double rate = 0.0;
	int heIndex = 0, vIndex = 0;
	// Initialize the dispersion sum
	int nHeSquare = 0, nVSquare = 0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Create the key for the maps
			auto key = std::make_pair(heIndex, vIndex);

			// Check if this cluster exists
			if (reactingMap.find(key) == reactingMap.end())
				continue;

			// Compute nSquare for the dispersion
			nHeSquare += heIndex * heIndex;
			nVSquare += vIndex * vIndex;

			// Get all the reaction vectors at this index
			reactingPairs = reactingMap[key];
			combiningReactants = combiningMap[key];
			dissociatingPairs = dissociatingMap[key];
			emissionPairs = emissionMap[key];

			// Initialize all the effective vectors
			effReactingPairs.clear();
			effCombiningReactants.clear();
			effDissociatingPairs.clear();
			effEmissionPairs.clear();

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
				reactingMap[key][i].kConstant = rate
						/ (double) nTot;

				// Add the reacting pair to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effReactingPairs.push_back(&reactingMap[key][i]);
				}
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
				combiningMap[key][i].kConstant = rate
						/ (double) nTot;

				// Add the combining reactant to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effCombiningReactants.push_back(
							&combiningMap[key][i]);

					// Add itself to the list again to account for the correct rate
					if (id == combiningReactant->getId())
						effCombiningReactants.push_back(
								&combiningMap[key][i]);
				}
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
					rate = calculateDissociationConstant(*dissociatingCluster,
						*this, *otherEmittedCluster);
				} else {
					// otherEmittedCluster is the single size one
					rate = calculateDissociationConstant(*dissociatingCluster,
							*otherEmittedCluster, *this);

				}
				// Set it in the pair
				dissociatingMap[key][i].kConstant = rate
					/ (double) nTot;

				// Add the dissociating pair to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effDissociatingPairs.push_back(
							&dissociatingMap[key][i]);

					// Add itself to the list again to account for the correct rate
					if (id == otherEmittedCluster->getId())
						effDissociatingPairs.push_back(
								&dissociatingMap[key][i]);
				}
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
				emissionMap[key][i].kConstant = rate
						/ (double) nTot;

				// Add the emission pair to the effective vector
				// if the rate is not 0.0
				if (!xolotlCore::equal(rate, 0.0)) {
					effEmissionPairs.push_back(&emissionMap[key][i]);
				}
			}

			// Shrink the arrays to save some space
			effReactingPairs.shrink_to_fit();
			effCombiningReactants.shrink_to_fit();
			effDissociatingPairs.shrink_to_fit();
			effEmissionPairs.shrink_to_fit();

			// Set the arrays in the effective maps
			effReactingMap[key] = effReactingPairs;
			effCombiningMap[key] = effCombiningReactants;
			effDissociatingMap[key] = effDissociatingPairs;
			effEmissionMap[key] = effEmissionPairs;
		}
	}

	// Reset the vectors to save memory
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Compute the dispersions
	dispersionHe = ((double) nHeSquare
			- (double) (compositionMap[heType] * compositionMap[heType])
					/ (double) nTot) / (double) nTot;
	dispersionV = ((double) nVSquare
			- (double) (compositionMap[vType] * compositionMap[vType])
					/ (double) nTot) / (double) nTot;

	if (equal(dispersionHe, 0.0))
		dispersionHe = 1.0;

	if (equal(dispersionV, 0.0))
		dispersionV = 1.0;

	return;
}

void SuperCluster::resetConnectivities() {
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

	// Local declarations
	int heIndex = 0, vIndex = 0;

	// Loop on the vacancy width
	for (int k = 0; k < sectionVWidth; k++) {
		// Compute the vacancy index
		vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

		// Loop on the helium width
		for (int j = 0; j < sectionHeWidth; j++) {
			// Compute the helium index
			heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

			// Create the key for the maps
			auto key = std::make_pair(heIndex, vIndex);

			// Check if this cluster exists
			if (effReactingMap.find(key) == effReactingMap.end())
				continue;

			// Get all the reaction vectors at this index
			effReactingPairs = effReactingMap[key];
			effCombiningReactants = effCombiningMap[key];
			effDissociatingPairs = effDissociatingMap[key];

			// Loop on the effective reacting pairs
			for (auto it = effReactingPairs.begin(); it != effReactingPairs.end();
					++it) {
				// The cluster is connecting to both clusters in the pair
				setReactionConnectivity((*it)->first->getId());
				setReactionConnectivity((*it)->first->getHeMomentumId());
				setReactionConnectivity((*it)->first->getVMomentumId());
				setReactionConnectivity((*it)->second->getId());
				setReactionConnectivity((*it)->second->getHeMomentumId());
				setReactionConnectivity((*it)->second->getVMomentumId());
			}

			// Loop on the effective combining reactants
			for (auto it = effCombiningReactants.begin();
					it != effCombiningReactants.end(); ++it) {
				// The cluster is connecting to the combining cluster
				setReactionConnectivity((*it)->combining->getId());
				setReactionConnectivity((*it)->combining->getHeMomentumId());
				setReactionConnectivity((*it)->combining->getVMomentumId());
			}

			// Loop on the effective dissociating pairs
			for (auto it = effDissociatingPairs.begin();
					it != effDissociatingPairs.end(); ++it) {
				// The cluster is connecting to the dissociating cluster which
				// is the first one by definition
				setDissociationConnectivity((*it)->first->getId());
				setDissociationConnectivity((*it)->first->getHeMomentumId());
				setDissociationConnectivity((*it)->first->getVMomentumId());
			}

			// Don't loop on the effective emission pairs because
			// this cluster is not connected to them
		}
	}

	return;
}

double SuperCluster::getDissociationFlux() {
	// Initial declarations
	int nPairs = 0, heIndex = 0, vIndex = 0;
	double flux = 0.0, heDistance = 0.0, vDistance = 0.0, heFactor = 0.0,
			vFactor = 0.0, value = 0.0;
	PSICluster *dissociatingCluster;

	// Loop on the effective map
	for (auto mapIt = effDissociatingMap.begin(); mapIt != effDissociatingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		heFactor = heDistance / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;
		vFactor = vDistance / dispersionV;

		// Get the pairs
		auto pairs = mapIt->second;
		// Set the total number of reactants that dissociate to form this one
		nPairs = pairs.size();
		// Loop over all dissociating clusters that form this cluster
		for (int i = 0; i < nPairs; i++) {
			// Get the dissociating cluster
			dissociatingCluster = pairs[i]->first;
			// Calculate the Dissociation flux
			value = pairs[i]->kConstant
					* dissociatingCluster->getConcentration(
							pairs[i]->firstHeDistance, pairs[i]->firstVDistance);
			flux += value;
			// Compute the momentum fluxes
			heMomentumFlux += value * heFactor;
			vMomentumFlux += value * vFactor;
		}
	}

	// Return the flux
	return flux;
}

double SuperCluster::getEmissionFlux() {
	// Initial declarations
	int nPairs = 0, heIndex = 0, vIndex = 0;
	double flux = 0.0, heDistance = 0.0, vDistance = 0.0, heFactor = 0.0,
			vFactor = 0.0, value = 0.0;

	// Loop on the effective map
	for (auto mapIt = effEmissionMap.begin(); mapIt != effEmissionMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		heFactor = heDistance / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;
		vFactor = vDistance / dispersionV;

		// Set the total number of reactants emitted by this one
		nPairs = mapIt->second.size();
		// Loop over all the pairs
		for (int i = 0; i < nPairs; i++) {
			// Update the flux
			value = mapIt->second[i]->kConstant * getConcentration(heDistance, vDistance);
			flux += value;
			// Compute the momentum fluxes
			heMomentumFlux -= value * heFactor;
			vMomentumFlux -= value * vFactor;
		}
	}

	return flux;
}

double SuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0, heDistance = 0.0, vDistance = 0.0, heFactor = 0.0,
			vFactor = 0.0, value = 0.0;
	int nPairs = 0, heIndex = 0, vIndex = 0;
	PSICluster *firstReactant, *secondReactant;

	// Loop on the effective map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		heFactor = heDistance / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;
		vFactor = vDistance / dispersionV;

		// Get the pairs
		auto pairs = mapIt->second;
		// Set the total number of reactants that produce to form this one
		nPairs = pairs.size();
		// Loop over all the reacting pairs
		for (int i = 0; i < nPairs; i++) {
			// Get the two reacting clusters
			firstReactant = pairs[i]->first;
			secondReactant = pairs[i]->second;
			// Update the flux
			value = pairs[i]->kConstant
					* firstReactant->getConcentration(pairs[i]->firstHeDistance, pairs[i]->firstVDistance)
					* secondReactant->getConcentration(pairs[i]->secondHeDistance, pairs[i]->secondVDistance);
			flux += value;
			// Compute the momentum fluxes
			heMomentumFlux += value * heFactor;
			vMomentumFlux += value * vFactor;
		}
	}

	// Return the production flux
	return flux;
}

double SuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0, heDistance = 0.0, vDistance = 0.0, heFactor = 0.0,
			vFactor = 0.0, value = 0.0;
	int nReactants = 0, heIndex = 0, vIndex = 0;
	PSICluster *combiningCluster;

	// Loop on the effective map
	for (auto mapIt = effCombiningMap.begin(); mapIt != effCombiningMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		heFactor = heDistance / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;
		vFactor = vDistance / dispersionV;

		// Get the pairs
		auto reactants = mapIt->second;
		// Set the total number of reactants that combine to form this one
		nReactants = reactants.size();
		// Loop over all possible clusters
		for (int i = 0; i < nReactants; i++) {
			// Get the cluster that combines with this one
			combiningCluster = reactants[i]->combining;
			// Calculate Second term of production flux
			value = reactants[i]->kConstant
					* combiningCluster->getConcentration(reactants[i]->heDistance, reactants[i]->vDistance)
					* getConcentration(heDistance, vDistance);
			flux += value;
			// Compute the momentum fluxes
			heMomentumFlux -= value * heFactor;
			vMomentumFlux -= value * vFactor;
		}
	}

	return flux;
}

void SuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, index = 0, heIndex = 0, vIndex = 0;
	double value = 0.0;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop on the effective map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		// Compute the vacancy index
		vIndex = mapIt->first.second;

		// Get the pairs
		auto pairs = mapIt->second;
		numReactants = pairs.size();
		for (int i = 0; i < numReactants; i++) {
			// Compute the contribution from the first part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->second->getConcentration(pairs[i]->secondHeDistance,
							pairs[i]->secondVDistance);
			index = pairs[i]->first->getId() - 1;
			partials[index] += value;
			index = pairs[i]->first->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->firstHeDistance;
			index = pairs[i]->first->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->firstVDistance;
			// Compute the contribution from the second part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->first->getConcentration(pairs[i]->firstHeDistance,
							pairs[i]->firstVDistance);
			index = pairs[i]->second->getId() - 1;
			partials[index] += value;
			index = pairs[i]->second->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->secondHeDistance;
			index = pairs[i]->second->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->secondVDistance;
		}
	}

	return;
}

void SuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, otherIndex = 0, heIndex = 0, vIndex = 0;
	PSICluster *cluster;
	double value = 0.0, heDistance = 0.0, vDistance = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Loop on the effective map
	for (auto mapIt = effCombiningMap.begin(); mapIt != effCombiningMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;

		// Get the pairs
		auto reactants = mapIt->second;
		numReactants = reactants.size();
		for (int i = 0; i < numReactants; i++) {
			cluster = (PSICluster *) reactants[i]->combining;
			// Remember that the flux due to combinations is OUTGOING (-=)!
			// Compute the contribution from this cluster
			value = reactants[i]->kConstant
					* cluster->getConcentration(reactants[i]->heDistance,
							reactants[i]->vDistance);
			partials[id - 1] -= value;
			partials[heMomId - 1] -= value * heDistance;
			partials[vMomId - 1] -= value * vDistance;
			// Compute the contribution from the combining cluster
			value = reactants[i]->kConstant
					* getConcentration(heDistance, vDistance);
			otherIndex = cluster->getId() - 1;
			partials[otherIndex] -= value;
			otherIndex = cluster->getHeMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->heDistance;
			otherIndex = cluster->getVMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->vDistance;
		}
	}

	return;
}

void SuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0, heIndex = 0, vIndex = 0;
	PSICluster *cluster;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Loop on the effective map
	for (auto mapIt = effDissociatingMap.begin(); mapIt != effDissociatingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		// Compute the vacancy index
		vIndex = mapIt->first.second;

		// Get the pairs
		auto pairs = mapIt->second;
		numPairs = pairs.size();
		for (int i = 0; i < numPairs; i++) {
			// Get the dissociating cluster
			cluster = pairs[i]->first;
			value = pairs[i]->kConstant;
			index = cluster->getId() - 1;
			partials[index] += value;
			index = cluster->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->firstHeDistance;
			index = cluster->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->firstVDistance;
		}
	}

	return;
}

void SuperCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0, heIndex = 0, vIndex = 0;
	double value = 0.0, heDistance = 0.0, vDistance = 0.0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Loop on the effective map
	for (auto mapIt = effEmissionMap.begin(); mapIt != effEmissionMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;

		// Get the pairs
		auto pairs = mapIt->second;
		numPairs = pairs.size();
		for (int i = 0; i < numPairs; i++) {
			// Modify the partial derivative. Remember that the flux
			// due to emission is OUTGOING (-=)!
			value = pairs[i]->kConstant;
			partials[id - 1] -= value;
			partials[heMomId - 1] -= value * heDistance;
			partials[vMomId - 1] -= value * vDistance;
		}
	}

	return;
}

void SuperCluster::getHeMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// For the momentum, the partial derivatives are the same as before, multiplied by the
	// distance to the mean

	// Local declarations
	int nPairs = 0, index = 0, otherIndex = 0, heIndex = 0, vIndex = 0;
	PSICluster *cluster;
	double value = 0.0, heDistance = 0.0, vDistance = 0.0, heFactor = 0.0;

	// Loop on the effective map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		heFactor = heDistance / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;

		// Get the current key
		auto key = mapIt->first;

		// Get the pairs
		auto pairs = mapIt->second;
		// Set the total number of reactants producing this one
		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Compute the contribution from the first part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->second->getConcentration(pairs[i]->secondHeDistance,
							pairs[i]->secondVDistance) * heFactor;
			index = pairs[i]->first->getId() - 1;
			partials[index] += value;
			index = pairs[i]->first->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->firstHeDistance;
			index = pairs[i]->first->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->firstVDistance;
			// Compute the contribution from the second part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->first->getConcentration(pairs[i]->firstHeDistance,
							pairs[i]->firstVDistance)
					* heFactor;
			index = pairs[i]->second->getId() - 1;
			partials[index] += value;
			index = pairs[i]->second->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->secondHeDistance;
			index = pairs[i]->second->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->secondVDistance;
		}

		// Get all the effective combining reactants at this index
		auto reactants = effCombiningMap.at(key);

		nPairs = reactants.size();
		for (int i = 0; i < nPairs; i++) {
			cluster = (PSICluster *) reactants[i]->combining;
			// Remember that the flux due to combinations is OUTGOING (-=)!
			// Compute the contribution from this cluster
			value = reactants[i]->kConstant
					* cluster->getConcentration(reactants[i]->heDistance,
							reactants[i]->vDistance)
					* heFactor;
			partials[id - 1] -= value;
			partials[heMomId - 1] -= value * heDistance;
			partials[vMomId - 1] -= value * vDistance;
			// Compute the contribution from the combining cluster
			value = reactants[i]->kConstant * getConcentration(heDistance, vDistance)
					* heFactor;
			otherIndex = cluster->getId() - 1;
			partials[otherIndex] -= value;
			otherIndex = cluster->getHeMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->heDistance;
			otherIndex = cluster->getVMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->vDistance;
		}

		// Get all the effective dissociating pairs at this index
		pairs = effDissociatingMap.at(key);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Get the dissociating cluster
			cluster = pairs[i]->first;
			value = pairs[i]->kConstant * heFactor;
			index = cluster->getId() - 1;
			partials[index] += value;
			index = cluster->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->firstHeDistance;
			index = cluster->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->firstVDistance;
		}

		// Get all the effective emission pairs at this index
		pairs = effEmissionMap.at(key);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Modify the partial derivative. Remember that the flux
			// due to emission is OUTGOING (-=)!
			value = pairs[i]->kConstant * heFactor;
			partials[id - 1] -= value;
			partials[heMomId - 1] -= value * heDistance;
			partials[vMomId - 1] -= value * vDistance;
		}
	}

	return;
}

void SuperCluster::getVMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// For the momentum, the partial derivatives are the same as before, multiplied by the
	// distance to the mean

	// Local declarations
	int nPairs = 0, index = 0, otherIndex = 0, heIndex = 0, vIndex = 0;
	PSICluster *cluster;
	double value = 0.0, heDistance = 0.0, vDistance = 0.0, vFactor = 0.0;

	// Loop on the effective map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = (double) heIndex - numHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = (double) vIndex - numV;
		vFactor = vDistance / dispersionV;

		// Get the current key
		auto key = mapIt->first;

		// Get the pairs
		auto pairs = mapIt->second;
		// Set the total number of reactants producing this one
		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Compute the contribution from the first part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->second->getConcentration(pairs[i]->secondHeDistance,
							pairs[i]->secondVDistance) * vFactor;
			index = pairs[i]->first->getId() - 1;
			partials[index] += value;
			index = pairs[i]->first->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->firstHeDistance;
			index = pairs[i]->first->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->firstVDistance;
			// Compute the contribution from the second part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->first->getConcentration(pairs[i]->firstHeDistance,
							pairs[i]->firstVDistance)
					* vFactor;
			index = pairs[i]->second->getId() - 1;
			partials[index] += value;
			index = pairs[i]->second->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->secondHeDistance;
			index = pairs[i]->second->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->secondVDistance;
		}

		// Get all the effective combining reactants at this index
		auto reactants = effCombiningMap.at(key);

		nPairs = reactants.size();
		for (int i = 0; i < nPairs; i++) {
			cluster = (PSICluster *) reactants[i]->combining;
			// Remember that the flux due to combinations is OUTGOING (-=)!
			// Compute the contribution from this cluster
			value = reactants[i]->kConstant
					* cluster->getConcentration(reactants[i]->heDistance,
							reactants[i]->vDistance)
					* vFactor;
			partials[id - 1] -= value;
			partials[heMomId - 1] -= value * heDistance;
			partials[vMomId - 1] -= value * vDistance;
			// Compute the contribution from the combining cluster
			value = reactants[i]->kConstant * getConcentration(heDistance, vDistance)
					* vFactor;
			otherIndex = cluster->getId() - 1;
			partials[otherIndex] -= value;
			otherIndex = cluster->getHeMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->heDistance;
			otherIndex = cluster->getVMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->vDistance;
		}

		// Get all the effective dissociating pairs at this index
		pairs = effDissociatingMap.at(key);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Get the dissociating cluster
			cluster = pairs[i]->first;
			value = pairs[i]->kConstant * vFactor;
			index = cluster->getId() - 1;
			partials[index] += value;
			index = cluster->getHeMomentumId() - 1;
			partials[index] += value * pairs[i]->firstHeDistance;
			index = cluster->getVMomentumId() - 1;
			partials[index] += value * pairs[i]->firstVDistance;
		}

		// Get all the effective emission pairs at this index
		pairs = effEmissionMap.at(key);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Modify the partial derivative. Remember that the flux
			// due to emission is OUTGOING (-=)!
			value = pairs[i]->kConstant * vFactor;
			partials[id - 1] -= value;
			partials[heMomId - 1] -= value * heDistance;
			partials[vMomId - 1] -= value * vDistance;
		}
	}

	return;
}
