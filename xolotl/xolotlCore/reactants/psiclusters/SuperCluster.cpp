// Includes
#include "SuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

/**
 * The helium momentum partials. It is computed only for super clusters.
 */
std::vector<double> heMomentumPartials;

/**
 * The vacancy momentum partials. It is computed only for super clusters.
 */
std::vector<double> vMomentumPartials;

SuperCluster::SuperCluster(double numHe, double numV, int nTot, int heWidth,
		int vWidth, double radius, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry), numHe(numHe), numV(numV), nTot(nTot), l0(0.0), l1He(
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

SuperCluster::SuperCluster(SuperCluster &other) :
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
	effReactingList = other.effReactingList;
	effCombiningList = other.effCombiningList;
	effDissociatingList = other.effDissociatingList;
	effEmissionList = other.effEmissionList;

	return;
}

std::shared_ptr<IReactant> SuperCluster::clone() {
	std::shared_ptr<IReactant> reactant(new SuperCluster(*this));

	return reactant;
}

double SuperCluster::getConcentration(double distHe, double distV) const {
	return l0 + (distHe * l1He) + (distV * l1V);
}

double SuperCluster::getHeMomentum() const {
	return l1He;
}

double SuperCluster::getVMomentum() const {
	return l1V;
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
			heDistance = 2.0 * (double) (heIndex - numHe) / (double) sectionHeWidth;
			vDistance = 2.0 * (double) (vIndex - numV) / (double) sectionVWidth;

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
			heDistance = 2.0 * (double) (heIndex - numHe) / (double) sectionHeWidth;
			vDistance = 2.0 * (double) (vIndex - numV) / (double) sectionVWidth;

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) heIndex;
		}
	}

	return conc;
}

double SuperCluster::getHeDistance(int he) const {
	return 2.0 * (double) (he - numHe) / (double) sectionHeWidth;
}

double SuperCluster::getVDistance(int v) const {
	return 2.0 * (double) (v - numV) / (double) sectionVWidth;
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
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;
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

					// Check if the rate is the biggest one up to now
					if (rate > biggestProductionRate)
						biggestProductionRate = rate;
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
				// otherEmittedCluster is the single size one
				rate = calculateDissociationConstant(*dissociatingCluster,
						*otherEmittedCluster, *this);

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

	// Set the biggest rate to the biggest production rate
	biggestRate = biggestProductionRate;

	// Compute the dispersions
	dispersionHe = 2.0 * ((double) nHeSquare
			- (double) (compositionMap[heType] * compositionMap[heType])
					/ (double) nTot) / ((double) (nTot * sectionHeWidth));
	dispersionV = 2.0 * ((double) nVSquare
			- (double) (compositionMap[vType] * compositionMap[vType])
					/ (double) nTot) / ((double) (nTot * sectionVWidth));

	if (equal(dispersionHe, 0.0))
		dispersionHe = 1.0;

	if (equal(dispersionV, 0.0))
		dispersionV = 1.0;

	// Method to optimize the reaction vectors
	optimizeReactions();

	return;
}

void SuperCluster::optimizeReactions() {
	// Local declarations
	double heFactor = 0.0, vFactor = 0.0, heDistance = 0.0, vDistance = 0.0;
	int nPairs = 0;
	PSICluster *firstReactant, *secondReactant, *combiningReactant,
			*dissociatingCluster, *otherEmittedCluster, *firstCluster,
			*secondCluster;
	int heIndex = 0, vIndex = 0;

	// Loop on the effective reacting map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end(); ) {
			// Get the two reacting clusters
			firstReactant = (*it)->first;
			secondReactant = (*it)->second;

			// Create a new SuperClusterProductionPair
			SuperClusterProductionPair superPair(firstReactant, secondReactant, (*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effReactingMap.end(); ++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end(); ) {
					// Get the two reacting clusters
					auto firstReactantBis = (*itBis)->first;
					auto secondReactantBis = (*itBis)->second;

					// Check if it is the same reaction
					if (firstReactantBis == firstReactant && secondReactantBis == secondReactant) {
						superPair.a00000 += 1.0;
						superPair.a00001 += heFactor;
						superPair.a00010 += vFactor;
						superPair.a00011 += (*itBis)->firstHeDistance;
						superPair.a00100 += (*itBis)->firstHeDistance * heFactor;
						superPair.a00101 += (*itBis)->firstHeDistance * vFactor;
						superPair.a00110 += (*itBis)->firstVDistance;
						superPair.a00111 += (*itBis)->firstVDistance * heFactor;
						superPair.a01000 += (*itBis)->firstVDistance * vFactor;
						superPair.a01001 += (*itBis)->secondHeDistance;
						superPair.a01010 += (*itBis)->secondHeDistance * heFactor;
						superPair.a01011 += (*itBis)->secondHeDistance * vFactor;
						superPair.a01100 += (*itBis)->secondVDistance;
						superPair.a01101 += (*itBis)->secondVDistance * heFactor;
						superPair.a01110 += (*itBis)->secondVDistance * vFactor;
						superPair.a01111 += (*itBis)->firstHeDistance * (*itBis)->secondHeDistance;
						superPair.a10000 += (*itBis)->firstHeDistance * (*itBis)->secondHeDistance * heFactor;
						superPair.a10001 += (*itBis)->firstHeDistance * (*itBis)->secondHeDistance * vFactor;
						superPair.a10010 += (*itBis)->firstHeDistance * (*itBis)->secondVDistance;
						superPair.a10011 += (*itBis)->firstHeDistance * (*itBis)->secondVDistance * heFactor;
						superPair.a10100 += (*itBis)->firstHeDistance * (*itBis)->secondVDistance * vFactor;
						superPair.a10101 += (*itBis)->firstVDistance * (*itBis)->secondHeDistance;
						superPair.a10110 += (*itBis)->firstVDistance * (*itBis)->secondHeDistance * heFactor;
						superPair.a10111 += (*itBis)->firstVDistance * (*itBis)->secondHeDistance * vFactor;
						superPair.a11000 += (*itBis)->firstVDistance * (*itBis)->secondVDistance;
						superPair.a11001 += (*itBis)->firstVDistance * (*itBis)->secondVDistance * heFactor;
						superPair.a11010 += (*itBis)->firstVDistance * (*itBis)->secondVDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else ++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effReactingList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Loop on the effective combining map
	for (auto mapIt = effCombiningMap.begin(); mapIt != effCombiningMap.end(); ++mapIt) {
		// Get the pairs
		auto clusters = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = clusters.begin(); it != clusters.end(); ) {
			// Get the combining cluster
			combiningReactant = (*it)->combining;

			// Create a new SuperClusterProductionPair with NULL as the second cluster because
			// we do not need it
			SuperClusterProductionPair superPair(combiningReactant, NULL, (*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effCombiningMap.end(); ++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heDistance = 2.0 * (double) (heIndex - numHe) / (double) sectionHeWidth;
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vDistance = 2.0 * (double) (vIndex - numV) / (double) sectionVWidth;
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto clustersBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = clustersBis.begin(); itBis != clustersBis.end(); ) {
					// Get the two reacting clusters
					auto combiningReactantBis = (*itBis)->combining;

					// Check if it is the same reaction
					if (combiningReactantBis == combiningReactant) {
						superPair.a00000 += 1.0;
						superPair.a00001 += heFactor;
						superPair.a00010 += vFactor;
						superPair.a00011 += (*itBis)->heDistance;
						superPair.a00100 += (*itBis)->heDistance * heFactor;
						superPair.a00101 += (*itBis)->heDistance * vFactor;
						superPair.a00110 += (*itBis)->vDistance;
						superPair.a00111 += (*itBis)->vDistance * heFactor;
						superPair.a01000 += (*itBis)->vDistance * vFactor;
						superPair.a01001 += heDistance;
						superPair.a01010 += heDistance * heFactor;
						superPair.a01011 += heDistance * vFactor;
						superPair.a01100 += vDistance;
						superPair.a01101 += vDistance * heFactor;
						superPair.a01110 += vDistance * vFactor;
						superPair.a01111 += (*itBis)->heDistance * heDistance;
						superPair.a10000 += (*itBis)->heDistance * heDistance * heFactor;
						superPair.a10001 += (*itBis)->heDistance * heDistance * vFactor;
						superPair.a10010 += (*itBis)->heDistance * vDistance;
						superPair.a10011 += (*itBis)->heDistance * vDistance * heFactor;
						superPair.a10100 += (*itBis)->heDistance * vDistance * vFactor;
						superPair.a10101 += (*itBis)->vDistance * heDistance;
						superPair.a10110 += (*itBis)->vDistance * heDistance * heFactor;
						superPair.a10111 += (*itBis)->vDistance * heDistance * vFactor;
						superPair.a11000 += (*itBis)->vDistance * vDistance;
						superPair.a11001 += (*itBis)->vDistance * vDistance * heFactor;
						superPair.a11010 += (*itBis)->vDistance * vDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = clustersBis.erase(itBis);
					}
					// Go to the next element
					else ++itBis;
				}

				// Give back the pairs
				mapItBis->second = clustersBis;
			}

			// Add the super pair
			effCombiningList.push_front(superPair);

			// Remove the reaction from the vector
			it = clusters.erase(it);
		}
	}

	// Loop on the effective dissociating map
	for (auto mapIt = effDissociatingMap.begin(); mapIt != effDissociatingMap.end(); ++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end(); ) {
			// Get the two reacting clusters
			dissociatingCluster = (*it)->first;
			otherEmittedCluster = (*it)->second;

			// Create a new SuperClusterProductionPair
			SuperClusterDissociationPair superPair(dissociatingCluster, otherEmittedCluster, (*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effDissociatingMap.end(); ++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end(); ) {
					// Get the two reacting clusters
					auto dissociatingClusterBis = (*itBis)->first;
					auto otherEmittedClusterBis = (*itBis)->second;

					// Check if it is the same reaction
					if (dissociatingClusterBis == dissociatingCluster && otherEmittedClusterBis == otherEmittedCluster) {
						superPair.a0000 += 1.0;
						superPair.a0001 += heFactor;
						superPair.a0010 += vFactor;
						superPair.a0011 += (*itBis)->firstHeDistance;
						superPair.a0100 += (*itBis)->firstHeDistance * heFactor;
						superPair.a0101 += (*itBis)->firstHeDistance * vFactor;
						superPair.a0110 += (*itBis)->firstVDistance;
						superPair.a0111 += (*itBis)->firstVDistance * heFactor;
						superPair.a1000 += (*itBis)->firstVDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else ++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effDissociatingList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Loop on the effective emission map
	for (auto mapIt = effEmissionMap.begin(); mapIt != effEmissionMap.end(); ++mapIt) {
		// Get the pairs
		auto pairs = mapIt->second;
		// Loop over all the reacting pairs
		for (auto it = pairs.begin(); it != pairs.end(); ) {
			// Get the two reacting clusters
			firstCluster = (*it)->first;
			secondCluster = (*it)->second;

			// Create a new SuperClusterProductionPair
			SuperClusterDissociationPair superPair(firstCluster, secondCluster, (*it)->kConstant);

			// Loop on the whole super cluster to fill this super pair
			for (auto mapItBis = mapIt; mapItBis != effEmissionMap.end(); ++mapItBis) {
				// Compute the helium index
				heIndex = mapItBis->first.first;
				heDistance = 2.0 * (double) (heIndex - numHe) / (double) sectionHeWidth;
				heFactor = (double) (heIndex - numHe) / dispersionHe;
				// Compute the vacancy index
				vIndex = mapItBis->first.second;
				vDistance = 2.0 * (double) (vIndex - numV) / (double) sectionVWidth;
				vFactor = (double) (vIndex - numV) / dispersionV;

				// Get the pairs
				auto pairsBis = mapItBis->second;
				// Set the total number of reactants that produce to form this one
				// Loop over all the reacting pairs
				for (auto itBis = pairsBis.begin(); itBis != pairsBis.end(); ) {
					// Get the two reacting clusters
					auto firstClusterBis = (*itBis)->first;
					auto secondClusterBis = (*itBis)->second;

					// Check if it is the same reaction
					if (firstClusterBis == firstCluster && secondClusterBis == secondCluster) {
						superPair.a0000 += 1.0;
						superPair.a0001 += heFactor;
						superPair.a0010 += vFactor;
						superPair.a0011 += heDistance;
						superPair.a0100 += heDistance * heFactor;
						superPair.a0101 += heDistance * vFactor;
						superPair.a0110 += vDistance;
						superPair.a0111 += vDistance * heFactor;
						superPair.a1000 += vDistance * vFactor;

						// Do not delete the element if it is the original one
						if (itBis == it) {
							++itBis;
							continue;
						}

						// Remove the reaction from the vector
						itBis = pairsBis.erase(itBis);
					}
					// Go to the next element
					else ++itBis;
				}

				// Give back the pairs
				mapItBis->second = pairsBis;
			}

			// Add the super pair
			effEmissionList.push_front(superPair);

			// Remove the reaction from the vector
			it = pairs.erase(it);
		}
	}

	// Clear the maps because they won't be used anymore
	effReactingPairs.clear();
	effCombiningReactants.clear();
	effDissociatingPairs.clear();
	effEmissionPairs.clear();
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

void SuperCluster::updateRateConstants() {
	// Local declarations
	PSICluster *firstReactant, *secondReactant, *combiningReactant,
		*dissociatingCluster, *otherEmittedCluster, *firstCluster,
		*secondCluster;
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on the reacting list
	for (auto it = effReactingList.begin(); it != effReactingList.end(); ++it) {
		// Get the reactants
		firstReactant = (*it).first;
		secondReactant = (*it).second;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*firstReactant,
				*secondReactant);
		// Set it in the pair
		(*it).kConstant = rate / (double) nTot;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on the combining list
	for (auto it = effCombiningList.begin(); it != effCombiningList.end(); ++it) {
		// Get the reactants
		combiningReactant = (*it).first;
		// Compute the reaction constant
		rate = calculateReactionRateConstant(*this, *combiningReactant);
		// Set it in the combining cluster
		(*it).kConstant = rate / (double) nTot;
	}

	// Loop on the dissociating list
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end(); ++it) {
		dissociatingCluster = (*it).first;
		// The second element of the pair is the cluster that is also
		// emitted by the dissociation
		otherEmittedCluster = (*it).second;
		// Compute the dissociation constant
		// The order of the cluster is important here because of the binding
		// energy used in the computation. It is taken from the type of the first cluster
		// which must be the single one
		// otherEmittedCluster is the single size one
		rate = calculateDissociationConstant(*dissociatingCluster,
				*otherEmittedCluster, *this);
		// Set it in the pair
		(*it).kConstant = rate / (double) nTot;
	}

	// Loop on the emission list
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		firstCluster = (*it).first;
		secondCluster = (*it).second;
		// Compute the dissociation rate
		rate = calculateDissociationConstant(*this, *firstCluster,
				*secondCluster);
		// Set it in the pair
		(*it).kConstant = rate / (double) nTot;
	}

	// Set the biggest rate to the biggest production rate
	biggestRate = biggestProductionRate;

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
	for (auto it = effCombiningList.begin(); it != effCombiningList.end(); ++it) {
		// The cluster is connecting to the combining cluster
		setReactionConnectivity((*it).first->getId());
		setReactionConnectivity((*it).first->getHeMomentumId());
		setReactionConnectivity((*it).first->getVMomentumId());
	}

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end(); ++it) {
		// The cluster is connecting to the combining cluster
		setDissociationConnectivity((*it).first->getId());
		setDissociationConnectivity((*it).first->getHeMomentumId());
		setDissociationConnectivity((*it).first->getVMomentumId());
	}

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the momentum
	int dof = network->size() + 2 * network->getAll(superType).size();
	heMomentumPartials.resize(dof, 0.0);
	vMomentumPartials.resize(dof, 0.0);

	return;
}

double SuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;
	PSICluster *dissociatingCluster;

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end(); ++it) {
		// Get the dissociating clusters
		dissociatingCluster = (*it).first;
		double l0A = dissociatingCluster->getConcentration(0.0, 0.0);
		double lHeA = dissociatingCluster->getHeMomentum();
		double lVA = dissociatingCluster->getVMomentum();
		// Update the flux
		value = (*it).kConstant;
		flux += value * ((*it).a0000 * l0A
				+ (*it).a0011 * lHeA
				+ (*it).a0110 * lVA);
		// Compute the momentum fluxes
		heMomentumFlux += value * ((*it).a0001 * l0A
				+ (*it).a0100 * lHeA
				+ (*it).a0111 * lVA);
		vMomentumFlux += value * ((*it).a0010 * l0A
				+ (*it).a0101 * lHeA
				+ (*it).a1000 * lVA);
	}

	// Return the flux
	return flux;
}

double SuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0, value = 0.0;

	// Loop over all the emission pairs
	for (auto it = effEmissionList.begin(); it != effEmissionList.end(); ++it) {
		// Update the flux
		value = (*it).kConstant;
		flux += value * ((*it).a0000 * l0
				+ (*it).a0011 * l1He
				+ (*it).a0110 * l1V);
		// Compute the momentum fluxes
		heMomentumFlux -= value * ((*it).a0001 * l0
				+ (*it).a0100 * l1He
				+ (*it).a0111 * l1V);
		vMomentumFlux -= value * ((*it).a0010 * l0
				+ (*it).a0101 * l1He
				+ (*it).a1000 * l1V);
	}

	return flux;
}

double SuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0, value = 0.0;
	PSICluster *firstReactant, *secondReactant;

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
		value = (*it).kConstant;
		flux += value * ((*it).a00000 * l0A * l0B
				+ (*it).a01001 * l0A * lHeB
				+ (*it).a01100 * l0A * lVB
				+ (*it).a00011 * lHeA * l0B
				+ (*it).a01111 * lHeA * lHeB
				+ (*it).a10010 * lHeA * lVB
				+ (*it).a00110 * lVA * l0B
				+ (*it).a10101 * lVA * lHeB
				+ (*it).a11000 * lVA * lVB);
		// Compute the momentum fluxes
		heMomentumFlux += value * ((*it).a00001 * l0A * l0B
				+ (*it).a01010 * l0A * lHeB
				+ (*it).a01101 * l0A * lVB
				+ (*it).a00100 * lHeA * l0B
				+ (*it).a10000 * lHeA * lHeB
				+ (*it).a10011 * lHeA * lVB
				+ (*it).a00111 * lVA * l0B
				+ (*it).a10110 * lVA * lHeB
				+ (*it).a11001 * lVA * lVB);
		vMomentumFlux += value * ((*it).a00010 * l0A * l0B
				+ (*it).a01011 * l0A * lHeB
				+ (*it).a01110 * l0A * lVB
				+ (*it).a00101 * lHeA * l0B
				+ (*it).a10001 * lHeA * lHeB
				+ (*it).a10100 * lHeA * lVB
				+ (*it).a01000 * lVA * l0B
				+ (*it).a10111 * lVA * lHeB
				+ (*it).a11010 * lVA * lVB);
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

	// Loop over all the combining clusters
	for (auto it = effCombiningList.begin(); it != effCombiningList.end(); ++it) {
		// Get the two reacting clusters
		combiningCluster = (*it).first;
		double l0A = combiningCluster->getConcentration(0.0, 0.0);
		double lHeA = combiningCluster->getHeMomentum();
		double lVA = combiningCluster->getVMomentum();
		// Update the flux
		value = (*it).kConstant;
		flux += value * ((*it).a00000 * l0A * l0
				+ (*it).a01001 * l0A * l1He
				+ (*it).a01100 * l0A * l1V
				+ (*it).a00011 * lHeA * l0
				+ (*it).a01111 * lHeA * l1He
				+ (*it).a10010 * lHeA * l1V
				+ (*it).a00110 * lVA * l0
				+ (*it).a10101 * lVA * l1He
				+ (*it).a11000 * lVA * l1V);
		// Compute the momentum fluxes
		heMomentumFlux -= value * ((*it).a00001 * l0A * l0
				+ (*it).a01010 * l0A * l1He
				+ (*it).a01101 * l0A * l1V
				+ (*it).a00100 * lHeA * l0
				+ (*it).a10000 * lHeA * l1He
				+ (*it).a10011 * lHeA * l1V
				+ (*it).a00111 * lVA * l0
				+ (*it).a10110 * lVA * l1He
				+ (*it).a11001 * lVA * l1V);
		vMomentumFlux -= value * ((*it).a00010 * l0A * l0
				+ (*it).a01011 * l0A * l1He
				+ (*it).a01110 * l0A * l1V
				+ (*it).a00101 * lHeA * l0
				+ (*it).a10001 * lHeA * l1He
				+ (*it).a10100 * lHeA * l1V
				+ (*it).a01000 * lVA * l0
				+ (*it).a10111 * lVA * l1He
				+ (*it).a11010 * lVA * l1V);
	}

	return flux;
}

void SuperCluster::getPartialDerivatives(std::vector<double> & partials) const {
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

void SuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	double value = 0.0;
	int index = 0;
	PSICluster *firstReactant, *secondReactant;

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
		value = (*it).kConstant;
		index = firstReactant->getId() - 1;
		partials[index] += value * ((*it).a00000 * l0B
				+ (*it).a01001 * lHeB
				+ (*it).a01100 * lVB);
		heMomentumPartials[index] += value * ((*it).a00001 * l0B
				+ (*it).a01010 * lHeB
				+ (*it).a01101 * lVB);
		vMomentumPartials[index] += value * ((*it).a00010 * l0B
				+ (*it).a01011 * lHeB
				+ (*it).a01110 * lVB);
		index = firstReactant->getHeMomentumId() - 1;
		partials[index] += value * ((*it).a00011 * l0B
				+ (*it).a01111 * lHeB
				+ (*it).a10010 * lVB);
		heMomentumPartials[index] += value * ((*it).a00100 * l0B
				+ (*it).a10000 * lHeB
				+ (*it).a10011 * lVB);
		vMomentumPartials[index] += value * ((*it).a00101 * l0B
				+ (*it).a10001 * lHeB
				+ (*it).a10100 * lVB);
		index = firstReactant->getVMomentumId() - 1;
		partials[index] += value * ((*it).a00110 * l0B
				+ (*it).a10101 * lHeB
				+ (*it).a11000 * lVB);
		heMomentumPartials[index] += value * ((*it).a00111 * l0B
				+ (*it).a10110 * lHeB
				+ (*it).a11001 * lVB);
		vMomentumPartials[index] += value * ((*it).a01000 * l0B
				+ (*it).a10111 * lHeB
				+ (*it).a11010 * lVB);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant->getId() - 1;
		partials[index] += value * ((*it).a00000 * l0A
				+ (*it).a00011 * lHeA
				+ (*it).a00110 * lVA);
		heMomentumPartials[index] += value * ((*it).a00001 * l0A
				+ (*it).a00100 * lHeA
				+ (*it).a00111 * lVA);
		vMomentumPartials[index] += value * ((*it).a00010 * l0A
				+ (*it).a00101 * lHeA
				+ (*it).a01000 * lVA);
		index = secondReactant->getHeMomentumId() - 1;
		partials[index] += value * ((*it).a01001 * l0A
				+ (*it).a01111 * lHeA
				+ (*it).a10101 * lVA);
		heMomentumPartials[index] += value * ((*it).a01010 * l0A
				+ (*it).a10000 * lHeA
				+ (*it).a10110 * lVA);
		vMomentumPartials[index] += value * ((*it).a01011 * l0A
				+ (*it).a10001 * lHeA
				+ (*it).a10111 * lVA);
		index = secondReactant->getVMomentumId() - 1;
		partials[index] += value * ((*it).a01100 * l0A
				+ (*it).a10010 * lHeA
				+ (*it).a11000 * lVA);
		heMomentumPartials[index] += value * ((*it).a01101 * l0A
				+ (*it).a10011 * lHeA
				+ (*it).a11001 * lVA);
		vMomentumPartials[index] += value * ((*it).a01110 * l0A
				+ (*it).a10100 * lHeA
				+ (*it).a11010 * lVA);
	}

	return;
}

void SuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	PSICluster *cluster;
	double value = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Loop over all the combining clusters
		for (auto it = effCombiningList.begin(); it != effCombiningList.end(); ++it) {
		// Get the two reacting clusters
		cluster = (*it).first;
		double l0A = cluster->getConcentration(0.0, 0.0);
		double lHeA = cluster->getHeMomentum();
		double lVA = cluster->getVMomentum();

		// Compute the contribution from the combining cluster
		value = (*it).kConstant;
		index = cluster->getId() - 1;
		partials[index] -= value * ((*it).a00000 * l0
				+ (*it).a01001 * l1He
				+ (*it).a01100 * l1V);
		heMomentumPartials[index] -= value * ((*it).a00001 * l0
				+ (*it).a01010 * l1He
				+ (*it).a01101 * l1V);
		vMomentumPartials[index] -= value * ((*it).a00010 * l0
				+ (*it).a01011 * l1He
				+ (*it).a01110 * l1V);
		index = cluster->getHeMomentumId() - 1;
		partials[index] -= value * ((*it).a00011 * l0
				+ (*it).a01111 * l1He
				+ (*it).a10010 * l1V);
		heMomentumPartials[index] -= value * ((*it).a00100 * l0
				+ (*it).a10000 * l1He
				+ (*it).a10011 * l1V);
		vMomentumPartials[index] -= value * ((*it).a00101 * l0
				+ (*it).a10001 * l1He
				+ (*it).a10100 * l1V);
		index = cluster->getVMomentumId() - 1;
		partials[index] -= value * ((*it).a00110 * l0
				+ (*it).a10101 * l1He
				+ (*it).a11000 * l1V);
		heMomentumPartials[index] -= value * ((*it).a00111 * l0
				+ (*it).a10110 * l1He
				+ (*it).a11001 * l1V);
		vMomentumPartials[index] -= value * ((*it).a01000 * l0
				+ (*it).a10111 * l1He
				+ (*it).a11010 * l1V);
		// Compute the contribution from this cluster
		index = id - 1;
		partials[index] -= value * ((*it).a00000 * l0A
				+ (*it).a00011 * lHeA
				+ (*it).a00110 * lVA);
		heMomentumPartials[index] -= value * ((*it).a00001 * l0A
				+ (*it).a00100 * lHeA
				+ (*it).a00111 * lVA);
		vMomentumPartials[index] -= value * ((*it).a00010 * l0A
				+ (*it).a00101 * lHeA
				+ (*it).a01000 * lVA);
		index = heMomId - 1;
		partials[index] -= value * ((*it).a01001 * l0A
				+ (*it).a01111 * lHeA
				+ (*it).a10101 * lVA);
		heMomentumPartials[index] -= value * ((*it).a01010 * l0A
				+ (*it).a10000 * lHeA
				+ (*it).a10110 * lVA);
		vMomentumPartials[index] -= value * ((*it).a01011 * l0A
				+ (*it).a10001 * lHeA
				+ (*it).a10111 * lVA);
		index = vMomId - 1;
		partials[index] -= value * ((*it).a01100 * l0A
				+ (*it).a10010 * lHeA
				+ (*it).a11000 * lVA);
		heMomentumPartials[index] -= value * ((*it).a01101 * l0A
				+ (*it).a10011 * lHeA
				+ (*it).a11001 * lVA);
		vMomentumPartials[index] -= value * ((*it).a01110 * l0A
				+ (*it).a10100 * lHeA
				+ (*it).a11010 * lVA);
	}

	return;
}

void SuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int index = 0;
	PSICluster *cluster;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Loop over all the dissociating pairs
	for (auto it = effDissociatingList.begin(); it != effDissociatingList.end(); ++it) {
		// Get the dissociating clusters
		cluster = (*it).first;

		// Compute the contribution from the dissociating cluster
		value = (*it).kConstant;
		index = cluster->getId() - 1;
		partials[index] += value * ((*it).a0000);
		heMomentumPartials[index] += value * ((*it).a0001);
		vMomentumPartials[index] += value * ((*it).a0010);
		index = cluster->getHeMomentumId() - 1;
		partials[index] += value * ((*it).a0011);
		heMomentumPartials[index] += value * ((*it).a0100);
		vMomentumPartials[index] += value * ((*it).a0101);
		index = cluster->getVMomentumId() - 1;
		partials[index] += value * ((*it).a0110);
		heMomentumPartials[index] += value * ((*it).a0111);
		vMomentumPartials[index] += value * ((*it).a1000);
	}

	return;
}

void SuperCluster::getEmissionPartialDerivatives(
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
		value = (*it).kConstant;
		index = id - 1;
		partials[index] -= value * ((*it).a0000);
		heMomentumPartials[index] -= value * ((*it).a0001);
		vMomentumPartials[index] -= value * ((*it).a0010);
		index = heMomId - 1;
		partials[index] -= value * ((*it).a0011);
		heMomentumPartials[index] -= value * ((*it).a0100);
		vMomentumPartials[index] -= value * ((*it).a0101);
		index = vMomId - 1;
		partials[index] -= value * ((*it).a0110);
		heMomentumPartials[index] -= value * ((*it).a0111);
		vMomentumPartials[index] -= value * ((*it).a1000);
	}

	return;
}

void SuperCluster::getHeMomentPartialDerivatives(std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = heMomentumPartials[i];
	}

	return;
}

void SuperCluster::getVMomentPartialDerivatives(std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = vMomentumPartials[i];
	}

	return;
}

void SuperCluster::initializeBursting(int surfacePos,
		std::vector<double> grid) {
	// Add the needed reaction connectivity
	// Each V cluster connects to every HeV clusters with the same number of V

	// Initial declarations
	int heIndex = 0, vIndex = 0;

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

			// Get the corresponding V cluster to add the connectivity
			auto vCluster = (PSICluster *) network->get(vType, vIndex);
			vCluster->setReactionConnectivity(id);
			vCluster->setReactionConnectivity(heMomId);
			vCluster->setReactionConnectivity(vMomId);
		}
	}

	// Method that will fill the index vector that is actually used during the solving steps
	initializeBurstingIndex(surfacePos, grid);

	return;
}

void SuperCluster::initializeBurstingIndex(int surfacePos,
		std::vector<double> grid) {
	// Clear the vector of HeV bubble bursting at each grid point
	burstingIndexVector.clear();
	// Initial declarations
	int heIndex = 0, vIndex = 0;

	// Loop on the grid points
	for (int i = 0; i < grid.size(); i++) {
		// Boolean to know if the cluster bursts at this depth
		bool burst = false;

		// Get the depth
		double depth = grid[i] - grid[surfacePos];

		// Loop on the vacancy width
		for (int k = 0; k < sectionVWidth; k++) {
			// Compute the vacancy index
			vIndex = (int) (numV - (double) sectionVWidth / 2.0) + k + 1;

			// Loop on the helium width
			for (int j = 0; j < sectionHeWidth; j++) {
				// Compute the helium index
				heIndex = (int) (numHe - (double) sectionHeWidth / 2.0) + j + 1;

				auto pair = std::make_pair(heIndex, vIndex);

				// Check if this cluster exists
				if (effReactingMap.find(pair) == effReactingMap.end())
					continue;

				if (heIndex < vIndex * 3) continue;

				// Compute the radius
				double radius = (sqrt(3.0) / 4.0) * xolotlCore::latticeConstant
						+ pow((3.0 * pow(xolotlCore::latticeConstant, 3.0) * (double) vIndex)
										/ (8.0 * xolotlCore::pi), (1.0 / 3.0))
						- pow((3.0 * pow(xolotlCore::latticeConstant, 3.0))
										/ (8.0 * xolotlCore::pi), (1.0 / 3.0));

				// If the radius is bigger than the distance to the surface there is bursting
				if (radius > depth) {
					// Set burst to true
					burst = true;
				}
			}
		}

		// Add indices to the index vector
		burstingIndexVector.push_back(burst);
	}

	return;
}

void SuperCluster::computeBurstingFlux(int xi, double *updatedConcOffset,
		double kBursting) {
	// Initial declarations
	int heIndex = 0, vIndex = 0;
	double heDistance = 0.0, vDistance = 0.0, value = 0.0, heFactor = 0.0,
			vFactor = 0.0;

	// Check if this cluster bursts at this depth
	if (!burstingIndexVector[xi]) return;

	// Loop on the effective map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = 2.0 * (double) (heIndex - numHe) / (double) sectionHeWidth;
		heFactor = (double) (heIndex - numHe) / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = 2.0 * (double) (vIndex - numV) / (double) sectionVWidth;
		vFactor = (double) (vIndex - numV) / dispersionV;

		// Get the initial concentration
		double oldConc = getConcentration(heDistance, vDistance);

		// Get the V cluster with the same number of V
		auto vCluster = (PSICluster *) network->get(vType, vIndex);
		// And its ID
		int vId = vCluster->getId() - 1;

		// Update the concentrations (the bubble loses its concentration)
		value = kBursting * oldConc / (double) nTot;
		updatedConcOffset[id - 1] -= value;
		updatedConcOffset[heMomId - 1] -= value * heFactor;
		updatedConcOffset[vMomId - 1] -= value * vFactor;
		updatedConcOffset[vId] += value;
	}

	return;
}

int SuperCluster::computePartialsForBursting(double *val,
		int *indices, int xi, double kBursting, int iStart) {
	// Initial declarations
	int heIndex = 0, vIndex = 0, i = 0;
	double heDistance = 0.0, vDistance = 0.0, value = 0.0, heFactor = 0.0,
			vFactor = 0.0;

	// Check if this cluster bursts at this depth
	if (!burstingIndexVector[xi]) return 0;

	// Loop on the effective map
	for (auto mapIt = effReactingMap.begin(); mapIt != effReactingMap.end(); ++mapIt) {
		// Compute the helium index
		heIndex = mapIt->first.first;
		heDistance = 2.0 * (double) (heIndex - numHe) / (double) sectionHeWidth;
		heFactor = (double) (heIndex - numHe) / dispersionHe;
		// Compute the vacancy index
		vIndex = mapIt->first.second;
		vDistance = 2.0 * (double) (vIndex - numV) / (double) sectionVWidth;
		vFactor = (double) (vIndex - numV) / dispersionV;

		// Compute the value of the partial derivatives
		value = kBursting / (double) nTot;
		// Keep the super index
		indices[(iStart * 2) + (i * 24)] = id - 1;
		indices[(iStart * 2) + (i * 24) + 1] = id - 1;
		val[iStart + (i * 12)] = -value;
		indices[(iStart * 2) + (i * 24) + 2] = id - 1;
		indices[(iStart * 2) + (i * 24) + 3] = heMomId - 1;
		val[iStart + (i * 12) + 1] = -value * heDistance;
		indices[(iStart * 2) + (i * 24) + 4] = id - 1;
		indices[(iStart * 2) + (i * 24) + 5] = vMomId - 1;
		val[iStart + (i * 12) + 2] = -value * vDistance;

		// Set the He momentum partials
		indices[(iStart * 2) + (i * 24) + 6] = heMomId - 1;
		indices[(iStart * 2) + (i * 24) + 7] = id - 1;
		val[iStart + (i * 12) + 3] = -value * heFactor;
		indices[(iStart * 2) + (i * 24) + 8] = heMomId - 1;
		indices[(iStart * 2) + (i * 24) + 9] = heMomId - 1;
		val[iStart + (i * 12) + 4] = -value * heDistance * heFactor;
		indices[(iStart * 2) + (i * 24) + 10] = heMomId - 1;
		indices[(iStart * 2) + (i * 24) + 11] = vMomId - 1;
		val[iStart + (i * 12) + 5] = -value * vDistance * heFactor;

		// Set the V momentum partials
		indices[(iStart * 2) + (i * 24) + 12] = vMomId - 1;
		indices[(iStart * 2) + (i * 24) + 13] = id - 1;
		val[iStart + (i * 12) + 6] = -value * vFactor;
		indices[(iStart * 2) + (i * 24) + 14] = vMomId - 1;
		indices[(iStart * 2) + (i * 24) + 15] = heMomId - 1;
		val[iStart + (i * 12) + 7] = -value * heDistance * vFactor;
		indices[(iStart * 2) + (i * 24) + 16] = vMomId - 1;
		indices[(iStart * 2) + (i * 24) + 17] = vMomId - 1;
		val[iStart + (i * 12) + 8] = -value * vDistance * vFactor;

		// Get the V cluster with the same number of V
		auto vCluster = (PSICluster *) network->get(vType, vIndex);
		// And its ID
		int vId = vCluster->getId() - 1;
		// Set its partial derivative
		indices[(iStart * 2) + (i * 24) + 18] = vId;
		indices[(iStart * 2) + (i * 24) + 19] = id - 1;
		val[iStart + (i * 12) + 9] = value;
		indices[(iStart * 2) + (i * 24) + 20] = vId;
		indices[(iStart * 2) + (i * 24) + 21] = heMomId - 1;
		val[iStart + (i * 12) + 10] = value * heDistance;
		indices[(iStart * 2) + (i * 24) + 22] = vId;
		indices[(iStart * 2) + (i * 24) + 23] = vMomId - 1;
		val[iStart + (i * 12) + 11] = value * vDistance;

		// Increment i
		i++;
	}

	return i;
}

int SuperCluster::getNBursting(int xi) {
	return nTot * (int) burstingIndexVector[xi];
}
