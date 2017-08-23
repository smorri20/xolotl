#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <cassert>

using namespace xolotlCore;


ReactionNetwork::ReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		handlerRegistry(registry), temperature(0.0), networkSize(0), dissociationsEnabled(
				true), numVClusters(0), numIClusters(0), numSuperClusters(0), maxVClusterSize(
				0), maxIClusterSize(0) {
	// Counter for the number of times the network concentration is updated.
	concUpdateCounter = handlerRegistry->getEventCounter("net_conc_updates");
	// Setup the vector to hold all of the reactants

	return;
}


double ReactionNetwork::calculateReactionRateConstant(
		ProductionReaction * reaction) const {
	// Get the reaction radii
	double r_first = reaction->first->getReactionRadius();
	double r_second = reaction->second->getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = reaction->first->getDiffusionCoefficient();
	double secondDiffusion = reaction->second->getDiffusionCoefficient();

	// Calculate and return
	double k_plus = 4.0 * xolotlCore::pi
			* (r_first + r_second + xolotlCore::reactionRadius)
			* (firstDiffusion + secondDiffusion);
	return k_plus;
}

void ReactionNetwork::fillConcentrationsArray(double * concentrations) {

	// Fill the array
    for(auto& currReactant : getAll()) {
		auto id = currReactant->getId() - 1;
		concentrations[id] = currReactant->getConcentration();
	}

	return;
}

void ReactionNetwork::updateConcentrationsFromArray(double * concentrations) {
	// Local Declarations
	int id = -1;

	// Set the concentrations
	concUpdateCounter->increment();	// increment the update concentration counter
	for (auto iter = allReactants.begin(); iter != allReactants.end();
			++iter) {
		id = (*iter)->getId() - 1;
		(*iter)->setConcentration(concentrations[id]);
	}

	return;
}

void ReactionNetwork::setTemperature(double temp) {
	// Set the temperature
	temperature = temp;

	// Update the temperature for all of the clusters
	for (int i = 0; i < networkSize; i++) {
		// This part will set the temperature in each reactant
		// and recompute the diffusion coefficient
		allReactants.at(i)->setTemperature(temp);
	}

	return;
}

double ReactionNetwork::getTemperature() const {
	return temperature;
}

std::vector<IReactant *> ReactionNetwork::getAll(Species type) const {
	// Local Declarations
	std::vector<IReactant *> reactants;

	// Only pull the reactants if the name is valid
    // TODO fix clients of getAll so they can use directly without 
    // shared pointers.
    // TODO fix so can validate type string.
    reactants.reserve(clusterTypeMap.at(type).size());
    for (auto& currReactant : clusterTypeMap.at(type)) {

        reactants.push_back(currReactant.get());
    }

	return reactants;
}

std::shared_ptr<ProductionReaction> ReactionNetwork::addProductionReaction(
		std::shared_ptr<ProductionReaction> reaction) {
	// Check if the given ProductionReaction already exists.
	auto key = reaction->descriptiveKey();
	auto iter = productionReactionMap.find(key);
	if (iter != productionReactionMap.end()) {
		// We already knew about the reaction, so return the one we
		// already had defined.
		return iter->second;
	}

	// We did not yet know about the given reaction.
	// Save it.
	productionReactionMap.emplace(key, reaction);
	allProductionReactions.emplace_back(reaction);

	return reaction;
}

std::shared_ptr<DissociationReaction> ReactionNetwork::addDissociationReaction(
		std::shared_ptr<DissociationReaction> reaction) {

	// Check if we already know about this reaction.
	auto key = reaction->descriptiveKey();
	auto iter = dissociationReactionMap.find(key);
	if (iter != dissociationReactionMap.end()) {
		// We already knew about the reaction.
		// Return the existing one.
		return iter->second;
	}

	// We did not yet know about the given reaction.
	// Add it, but also link it to its reverse reaction.
	// First, create the reverse reaction to get a pointer to it.
	auto reverseReaction = std::make_shared<ProductionReaction>(reaction->first,
			reaction->second);
	// Add this reverse reaction to our set of known reactions.
	reverseReaction = addProductionReaction(reverseReaction);

	// Indicate that the reverse reaction is the reverse reaction
	// to the newly-added dissociation reaction.
	reaction->reverseReaction = reverseReaction.get();

	// Add the dissociation reaction to our set of known reactions.
	dissociationReactionMap.emplace(key, reaction);
	allDissociationReactions.emplace_back(reaction);

	// Return the newly-added dissociation reaction.
	return reaction;
}

