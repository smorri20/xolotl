// Includes
#include "SuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

SuperCluster::SuperCluster(double numHe, double numV, int width, double radius,
		double energy,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(1, registry), numHe(numHe), numV(numV), dispersion(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	compositionMap[heType] = (int) (numHe * (double) width);
	compositionMap[vType] = (int) (numV * (double) width);

	// Set the width
	sectionWidth = width;

	// Set the reaction radius and formation energy
	reactionRadius = radius;
	formationEnergy = energy;

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

	return;
}

std::shared_ptr<Reactant> SuperCluster::clone() {
	std::shared_ptr<Reactant> reactant(new SuperCluster(*this));

	return reactant;
}

double SuperCluster::getConcentration(double id) const {
	return l0 + (id * l1);
}

double SuperCluster::getTotalConcentration() const {
	// Initial declarations
	int groupingIndex = 0;
	double groupingDistance = 0.0, conc = 0.0;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(groupingDistance) * (double) groupingIndex;
	}

	return conc;
}

double SuperCluster::getDistance(int he) const {
	return (double) he - numHe;
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

		// Set them in the super cluster map
		reactingMap[comp[heType]] = react;
		combiningMap[comp[heType]] = combi;
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

		// Set them in the super cluster map
		dissociatingMap[comp[heType]] = disso;
		emissionMap[comp[heType]] = emi;
	}

	return;
}

void SuperCluster::computeRateConstants() {
	// Local declarations
	PSICluster *firstReactant, *secondReactant, *combiningReactant,
		*dissociatingCluster, *otherEmittedCluster, *firstCluster,
		*secondCluster;
	double rate = 0.0;
	// Initialize the dispersion sum
	int nSquare = 0;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		int groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		// Compute nSquare for the dispersion
		nSquare += groupingIndex * groupingIndex;

		// Get all the reaction vectors at this index
		reactingPairs = reactingMap[groupingIndex];
		combiningReactants = combiningMap[groupingIndex];
		dissociatingPairs = dissociatingMap[groupingIndex];
		emissionPairs = emissionMap[groupingIndex];

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
			reactingMap[groupingIndex][i].kConstant = rate / (double) sectionWidth;

			// Add the reacting pair to the effective vector
			// if the rate is not 0.0
			if (!xolotlCore::equal(rate, 0.0)) {
				effReactingPairs.push_back(&reactingMap[groupingIndex][i]);
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
			combiningMap[groupingIndex][i].kConstant = rate / (double) sectionWidth;

			// Add the combining reactant to the effective vector
			// if the rate is not 0.0
			if (!xolotlCore::equal(rate, 0.0)) {
				effCombiningReactants.push_back(&combiningMap[groupingIndex][i]);

				// Add itself to the list again to account for the correct rate
				if (id == combiningReactant->getId())
					effCombiningReactants.push_back(&combiningMap[groupingIndex][i]);
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
				rate = calculateDissociationConstant(*dissociatingCluster, *this,
						*otherEmittedCluster);
			} else {
				// otherEmittedCluster is the single size one
				rate = calculateDissociationConstant(*dissociatingCluster,
						*otherEmittedCluster, *this);

			}
			// Set it in the pair
			dissociatingMap[groupingIndex][i].kConstant = rate / (double) sectionWidth;

			// Add the dissociating pair to the effective vector
			// if the rate is not 0.0
			if (!xolotlCore::equal(rate, 0.0)) {
				effDissociatingPairs.push_back(&dissociatingMap[groupingIndex][i]);

				// Add itself to the list again to account for the correct rate
				if (id == otherEmittedCluster->getId())
					effDissociatingPairs.push_back(&dissociatingMap[groupingIndex][i]);
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
			emissionMap[groupingIndex][i].kConstant = rate / (double) sectionWidth;

			// Add the emission pair to the effective vector
			// if the rate is not 0.0
			if (!xolotlCore::equal(rate, 0.0)) {
				effEmissionPairs.push_back(&emissionMap[groupingIndex][i]);
			}
		}

		// Shrink the arrays to save some space
		effReactingPairs.shrink_to_fit();
		effCombiningReactants.shrink_to_fit();
		effDissociatingPairs.shrink_to_fit();
		effEmissionPairs.shrink_to_fit();

		// Set the arrays in the effective maps
		effReactingMap[groupingIndex] = effReactingPairs;
		effCombiningMap[groupingIndex] = effCombiningReactants;
		effDissociatingMap[groupingIndex] = effDissociatingPairs;
		effEmissionMap[groupingIndex] = effEmissionPairs;
	}

	// Reset the vectors to save memory
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	// Compute the dispersion
	dispersion = ((double) nSquare - (double) (compositionMap[heType] * compositionMap[heType])
			/ (double) sectionWidth) / (double) sectionWidth;

	if (equal(dispersion, 0.0)) dispersion = 1.0;

	return;
}

void SuperCluster::resetConnectivities() {
	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(momId);
	setDissociationConnectivity(momId);

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		int groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		// Get all the reaction vectors at this index
		effReactingPairs = effReactingMap[groupingIndex];
		effCombiningReactants = effCombiningMap[groupingIndex];
		effDissociatingPairs = effDissociatingMap[groupingIndex];

		// Loop on the effective reacting pairs
		for (auto it = effReactingPairs.begin(); it != effReactingPairs.end(); ++it) {
			// The cluster is connecting to both clusters in the pair
			setReactionConnectivity((*it)->first->getId());
			setReactionConnectivity((*it)->first->getMomentumId());
			setReactionConnectivity((*it)->second->getId());
			setReactionConnectivity((*it)->second->getMomentumId());
		}

		// Loop on the effective combining reactants
		for (auto it = effCombiningReactants.begin(); it != effCombiningReactants.end(); ++it) {
			// The cluster is connecting to the combining cluster
			setReactionConnectivity((*it)->combining->getId());
			setReactionConnectivity((*it)->combining->getMomentumId());
		}

		// Loop on the effective dissociating pairs
		for (auto it = effDissociatingPairs.begin(); it != effDissociatingPairs.end(); ++it) {
			// The cluster is connecting to the dissociating cluster which
			// is the first one by definition
			setDissociationConnectivity((*it)->first->getId());
			setDissociationConnectivity((*it)->first->getMomentumId());
		}

		// Don't loop on the effective emission pairs because
		// this cluster is not connected to them
	}

	return;
}

double SuperCluster::getDissociationFlux() const {
	// Initial declarations
	int nPairs = 0, groupingIndex = 0;
	double flux = 0.0;
	PSICluster *dissociatingCluster;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		// Get all the effective dissociating pairs at this index
		auto pairs = effDissociatingMap.at(groupingIndex);

		// Set the total number of reactants that dissociate to form this one
		nPairs = pairs.size();
		// Loop over all dissociating clusters that form this cluster
		for (int i = 0; i < nPairs; i++) {
			// Get the dissociating cluster
			dissociatingCluster = pairs[i]->first;
			// Calculate the Dissociation flux
			flux += pairs[i]->kConstant
					* dissociatingCluster->getConcentration(pairs[i]->firstDistance);
		}
	}

	// Return the flux
	return flux;
}

double SuperCluster::getEmissionFlux() const {
	// Initial declarations
	int nPairs = 0, groupingIndex = 0;
	double flux = 0.0, groupingDistance = 0.0;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;
		// Get all the effective emission pairs at this index
		auto pairs = effEmissionMap.at(groupingIndex);

		// Set the total number of emission pairs
		nPairs = pairs.size();
		// Loop over all the pairs
		for (int i = 0; i < nPairs; i++) {
			// Update the flux
			flux += pairs[i]->kConstant * getConcentration(groupingDistance);
		}
	}

	return flux;
}

double SuperCluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	int nPairs = 0, groupingIndex = 0;
	PSICluster *firstReactant, *secondReactant;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		// Get all the effective reaction pairs at this index
		auto pairs = effReactingMap.at(groupingIndex);

		// Set the total number of reacting pairs
		nPairs = pairs.size();
		// Loop over all the reacting pairs
		for (int i = 0; i < nPairs; i++) {
			// Get the two reacting clusters
			firstReactant = pairs[i]->first;
			secondReactant = pairs[i]->second;
			// Update the flux
			flux += pairs[i]->kConstant * firstReactant->getConcentration(pairs[i]->firstDistance)
					* secondReactant->getConcentration(pairs[i]->secondDistance);
		}
	}

	// Return the production flux
	return flux;
}

double SuperCluster::getCombinationFlux() const {
	// Local declarations
	double flux = 0.0, groupingDistance = 0.0;
	int nReactants = 0, groupingIndex = 0;
	PSICluster *combiningCluster;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;
		// Get all the effective combining reactants at this index
		auto reactants = effCombiningMap.at(groupingIndex);

		// Set the total number of reactants that combine to form this one
		nReactants = reactants.size();
		// Loop over all possible clusters
		for (int i = 0; i < nReactants; i++) {
			// Get the cluster that combines with this one
			combiningCluster = reactants[i]->combining;
			// Calculate Second term of production flux
			flux += reactants[i]->kConstant
					* combiningCluster->getConcentration(reactants[i]->distance)
					* getConcentration(groupingDistance);
		}
	}

	return flux;
}

double SuperCluster::getMomentFlux() const {
	// For the momentum, the flux is the same as before, multiplied by the
	// distance to the mean

	// Local declarations
	double flux = 0.0, groupingDistance = 0.0, factor = 0.0;
	int nPairs = 0, groupingIndex = 0;
	PSICluster *firstReactant, *secondReactant, *dissociatingCluster, *combiningCluster;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index and the distance
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;
		factor = groupingDistance / dispersion;

		// Get all the effective dissociation pairs at this index
		auto pairs = effDissociatingMap.at(groupingIndex);

		// Set the total number of reactants that dissociate to form this one
		nPairs = pairs.size();
		// Loop over all dissociating clusters that form this cluster
		for (int i = 0; i < nPairs; i++) {
			// Get the dissociating cluster
			dissociatingCluster = pairs[i]->first;
			// Calculate the Dissociation flux
			flux += pairs[i]->kConstant
					* dissociatingCluster->getConcentration(pairs[i]->firstDistance) * factor;
		}

		// Get all the effective reacting pairs at this index
		pairs = effReactingMap.at(groupingIndex);

		// Set the total number of reacting pairs
		nPairs = pairs.size();
		// Loop over all the reacting pairs
		for (int i = 0; i < nPairs; i++) {
			// Get the two reacting clusters
			firstReactant = pairs[i]->first;
			secondReactant = pairs[i]->second;
			// Update the flux
			flux += pairs[i]->kConstant * firstReactant->getConcentration(pairs[i]->firstDistance)
					* secondReactant->getConcentration(pairs[i]->secondDistance) * factor;
		}

		// Get all the effective emission pairs at this index
		pairs = effEmissionMap.at(groupingIndex);

		// Set the total number of emission pairs
		nPairs = pairs.size();
		// Loop over all the pairs
		for (int i = 0; i < nPairs; i++) {
			// Update the flux
			flux -= pairs[i]->kConstant * getConcentration(groupingDistance) * factor;
		}

		// Get all the effective combining reactants at this index
		auto reactants = effCombiningMap.at(groupingIndex);

		// Set the total number of reactants that combine to form this one
		nPairs = reactants.size();
		// Loop over all possible clusters
		for (int i = 0; i < nPairs; i++) {
			// Get the cluster that combines with this one
			combiningCluster = reactants[i]->combining;
			// Calculate Second term of production flux
			flux -= reactants[i]->kConstant * combiningCluster->getConcentration(reactants[i]->distance)
					* getConcentration(groupingDistance) * factor;
		}
	}

	// Return the moment flux
	return flux;
}

void SuperCluster::getProductionPartialDerivatives(std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, index = 0, groupingIndex = 0;
	double value = 0.0;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		// Get all the effective reacting pairs at this index
		auto pairs = effReactingMap.at(groupingIndex);

		numReactants = pairs.size();
		for (int i = 0; i < numReactants; i++) {
			// Compute the contribution from the first part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->second->getConcentration(pairs[i]->secondDistance);
			index = pairs[i]->first->getId() - 1;
			partials[index] += value;
			index = pairs[i]->first->getMomentumId() - 1;
			partials[index] += value * pairs[i]->firstDistance;
			// Compute the contribution from the second part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->first->getConcentration(pairs[i]->firstDistance);
			index = pairs[i]->second->getId() - 1;
			partials[index] += value;
			index = pairs[i]->second->getMomentumId() - 1;
			partials[index] += value * pairs[i]->secondDistance;
		}
	}

	return;
}

void SuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, otherIndex = 0, groupingIndex = 0;
	PSICluster *cluster;
	double value = 0.0, groupingDistance = 0.0;

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;
		// Get all the effective combining reactants at this index
		auto reactants = effCombiningMap.at(groupingIndex);

		numReactants = reactants.size();
		for (int i = 0; i < numReactants; i++) {
			cluster = (PSICluster *) reactants[i]->combining;
			// Remember that the flux due to combinations is OUTGOING (-=)!
			// Compute the contribution from this cluster
			value = reactants[i]->kConstant
					* cluster->getConcentration(reactants[i]->distance);
			partials[id - 1] -= value;
			partials[momId - 1] -= value * groupingDistance;
			// Compute the contribution from the combining cluster
			value = reactants[i]->kConstant * getConcentration(groupingDistance);
			otherIndex = cluster->getId() - 1;
			partials[otherIndex] -= value;
			otherIndex = cluster->getMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->distance;
		}
	}

	return;
}

void SuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0, groupingIndex = 0;
	PSICluster *cluster;
	double value = 0.0;

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		// Get all the effective dissociating pairs at this index
		auto pairs = effDissociatingMap.at(groupingIndex);

		numPairs = pairs.size();
		for (int i = 0; i < numPairs; i++) {
			// Get the dissociating cluster
			cluster = pairs[i]->first;
			value = pairs[i]->kConstant;
			index = cluster->getId() - 1;
			partials[index] += value;
			index = cluster->getMomentumId() - 1;
			partials[index] += value * pairs[i]->firstDistance;
		}
	}

	return;
}

void SuperCluster::getEmissionPartialDerivatives(std::vector<double> & partials) const {
	// Initial declarations
	int numPairs = 0, index = 0, groupingIndex = 0;
	double value = 0.0, groupingDistance = 0.0;

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;
		// Get all the effective emission pairs at this index
		auto pairs = effEmissionMap.at(groupingIndex);

		numPairs = pairs.size();
		for (int i = 0; i < numPairs; i++) {
			// Modify the partial derivative. Remember that the flux
			// due to emission is OUTGOING (-=)!
			value = pairs[i]->kConstant;
			partials[id - 1] -= value;
			partials[momId - 1] -= value * groupingDistance;
		}
	}

	return;
}

void SuperCluster::getMomentPartialDerivatives(std::vector<double> & partials) const {
	// For the momentum, the partial derivatives are the same as before, multiplied by the
	// distance to the mean

	// Local declarations
	int nPairs = 0, index = 0, otherIndex = 0, groupingIndex = 0;
	PSICluster *cluster;
	double value = 0.0, groupingDistance = 0.0, factor = 0.0;

	// Loop on the width
	for (int j = 0; j < sectionWidth; j++) {
		// Compute the index and the distance
		groupingIndex = (int) (numHe - (double) sectionWidth / 2.0) + j + 1;
		groupingDistance = (double) groupingIndex - numHe;
		factor = groupingDistance / dispersion;

		// Get all the reaction vectors at this index
		auto pairs = effReactingMap.at(groupingIndex);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Compute the contribution from the first part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->second->getConcentration(pairs[i]->secondDistance) * factor;
			index = pairs[i]->first->getId() - 1;
			partials[index] += value;
			index = pairs[i]->first->getMomentumId() - 1;
			partials[index] += value * pairs[i]->firstDistance;
			// Compute the contribution from the second part of the reacting pair
			value = pairs[i]->kConstant
					* pairs[i]->first->getConcentration(pairs[i]->firstDistance) * factor;
			index = pairs[i]->second->getId() - 1;
			partials[index] += value;
			index = pairs[i]->second->getMomentumId() - 1;
			partials[index] += value * pairs[i]->secondDistance;
		}

		// Get all the effective combining reactants at this index
		auto reactants = effCombiningMap.at(groupingIndex);

		nPairs = reactants.size();
		for (int i = 0; i < nPairs; i++) {
			cluster = (PSICluster *) reactants[i]->combining;
			// Remember that the flux due to combinations is OUTGOING (-=)!
			// Compute the contribution from this cluster
			value = reactants[i]->kConstant
					* cluster->getConcentration(reactants[i]->distance) * factor;
			partials[id - 1] -= value;
			partials[momId - 1] -= value * groupingDistance;
			// Compute the contribution from the combining cluster
			value = reactants[i]->kConstant
					* getConcentration(groupingDistance) * factor;
			otherIndex = cluster->getId() - 1;
			partials[otherIndex] -= value;
			otherIndex = cluster->getMomentumId() - 1;
			partials[otherIndex] -= value * reactants[i]->distance;
		}

		// Get all the effective dissociating pairs at this index
		pairs = effDissociatingMap.at(groupingIndex);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Get the dissociating cluster
			cluster = pairs[i]->first;
			value = pairs[i]->kConstant * factor;
			index = cluster->getId() - 1;
			partials[index] += value;
			index = cluster->getMomentumId() - 1;
			partials[index] += value * pairs[i]->firstDistance;
		}

		// Get all the effective emission pairs at this index
		pairs = effEmissionMap.at(groupingIndex);

		nPairs = pairs.size();
		for (int i = 0; i < nPairs; i++) {
			// Modify the partial derivative. Remember that the flux
			// due to emission is OUTGOING (-=)!
			value = pairs[i]->kConstant * factor;
			partials[id - 1] -= value;
			partials[momId - 1] -= value * groupingDistance;
		}
	}

	return;
}
