#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <cassert>


namespace xolotlCore {

ReactionNetwork::ReactionNetwork(
        const std::set<ReactantType>& _knownReactantTypes,
        ReactantType _superClusterType,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> _registry) :
        knownReactantTypes(_knownReactantTypes),
        superClusterType(_superClusterType),
		handlerRegistry(_registry),
        temperature(0.0),
        dissociationsEnabled(true) {

	// Counter for the number of times the network concentration is updated.
	concUpdateCounter = handlerRegistry->getEventCounter("net_conc_updates");
	// Setup the vector to hold all of the reactants

    // Ensure our per-type cluster map can store Reactants of the types
    // we support.
    for (auto const& currType : knownReactantTypes) {
        clusterTypeMap.insert( {currType, IReactionNetwork::ReactantMap() } );
    }

    // Ensure we have a baseline for determining max cluster size for
    // the types we support.
    for (auto const& currType : knownReactantTypes) {
        maxClusterSizeMap.insert( {currType, 0 } );
    }
	return;
}


double ReactionNetwork::calculateReactionRateConstant(const ProductionReaction& reaction) const {

	// Get the reaction radii
	double r_first = reaction.first->getReactionRadius();
	double r_second = reaction.second->getReactionRadius();

	// Get the diffusion coefficients
	double firstDiffusion = reaction.first->getDiffusionCoefficient();
	double secondDiffusion = reaction.second->getDiffusionCoefficient();

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
	for (int i = 0; i < size(); i++) {
		// This part will set the temperature in each reactant
		// and recompute the diffusion coefficient
		allReactants.at(i)->setTemperature(temp);
	}

	return;
}

double ReactionNetwork::getTemperature() const {
	return temperature;
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

	// Return the newly-added dissociation reaction.
	return reaction;
}

void ReactionNetwork::dumpTo(std::ostream& os) const {

    // Dump flat view of reactants.
    os << size() << " reactants:\n";
    for (auto const& currReactant : allReactants) {
        os << *currReactant << '\n';
    }

    // Dump reactants of each type.
    // TODO what does this give us that the flat view doesn't?
    os << "per-type reactant map:\n";
    auto const& knownTypes = getKnownReactantTypes();
    for (auto const& currType : knownTypes) {
        auto const& currTypeReactantMap = clusterTypeMap.at(currType);
        os << currTypeReactantMap.size() << " " << toString(currType) << " reactants:\n";
        for (auto const& currMapItem : currTypeReactantMap) {
            os << *(currMapItem.second) << '\n';
        }
    }

    // Dump ProductionReactions.
    os << productionReactionMap.size() << " production reactions:\n";
    for (auto const& currMapItem : productionReactionMap) {
        os << *(currMapItem.second) << '\n';
    }

    // Dump DissociationReactions.
    os << dissociationReactionMap.size() << " dissociation reactions:\n";
    for (auto const& currMapItem : dissociationReactionMap) {
        os << *(currMapItem.second) << '\n';
    }
}

} // xolotlCore
