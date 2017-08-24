#include <sstream>
#include "ReactionNetwork.h"

namespace xolotlCore {

ProductionReaction::ProductionReaction(IReactant * _reactant1,
		IReactant * _reactant2) :
		first(_reactant1), second(_reactant2) {

	// We made an assumption about ordering.
	// Check if we were wrong, and if so, swap.
	if (_reactant2->getComposition() < _reactant1->getComposition()) {

        std::swap(first, second);
	}

    auto const& firstComp = first->getComposition();
    auto const& secondComp = second->getComposition();

	// Determine our descriptive key and cache it.
	// This assumes that our reactants are ordered.
    descKey.reserve(firstComp.size() + secondComp.size());
    descKey.insert(descKey.end(), firstComp.begin(), firstComp.end());
    descKey.insert(descKey.end(), secondComp.begin(), secondComp.end());
}

DissociationReaction::DissociationReaction(IReactant * dissociatingPtr,
		IReactant * _reactant1, IReactant * _reactant2) :
		dissociating(dissociatingPtr), first(_reactant1), second(_reactant2), reverseReaction(
				nullptr) {

	// We made an assumption about ordering.
	// Check if we were wrong, and if so, swap.
	if (_reactant2->getComposition() < _reactant1->getComposition()) {

        std::swap(first, second);
	}

    auto const& dissComp = dissociating->getComposition();
    auto const& firstComp = first->getComposition();
    auto const& secondComp = second->getComposition();

	// Determine our descriptive key and cache it.
	// This assumes that our reactants are ordered.
    descKey.reserve(dissComp.size() + firstComp.size() + secondComp.size());
    descKey.insert(descKey.end(), dissComp.begin(), dissComp.end());
    descKey.insert(descKey.end(), firstComp.begin(), firstComp.end());
    descKey.insert(descKey.end(), secondComp.begin(), secondComp.end());
}

} // namespace xolotlCore

