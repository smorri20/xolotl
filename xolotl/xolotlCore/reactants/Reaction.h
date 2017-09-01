#ifndef XCORE_REACTION_H
#define XCORE_REACTION_H

#include "xolotlCore/reactants/reactantsConfig.h"
#include "IReactant.h"

namespace xolotlCore {

#if defined(USE_ORIG_REACTANT_COMP_STRING)
// The original implementation ordered reactants first/second based
// on lexicographical ordering of their composition strings.
// We would prefer not to use composition strings for performance reasons,
// but for verification that we are initializing the same reaction network,
// we need to use the same comparison function.
inline
bool compStringCompare(const IReactant& _r1, const IReactant& _r2) {

    // build composition strings for the two reactant compositions,
    // that is in the same format as their old composition strings.
    std::string r1str = getCompString(_r1);
    std::string r2str = getCompString(_r2);

    return r1str < r2str;
}
#endif // defined(USE_ORIG_REACTANT_COMP_STRING)

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k is stored here and other information is implemented by
 * the daughter classes. k is computed when setTemperature() is called.
 */
class Reaction {
private:
    /**
     * Were the given reactant parameters correctly ordered when
     * we were constructed?
     * (Used only during construction to support detection of need to reorder.)
     */
    bool paramsCorrectlyOrdered;


protected:
    /**
     * Construct a Reaction.
     * @param _r1 One of the reactants.
     * @param _r2 The other reactant.
     */
    Reaction(IReactant* _r1, IReactant* _r2)
#if defined(USE_ORIG_REACTANT_COMP_STRING)
      : paramsCorrectlyOrdered( compStringCompare(*_r1, *_r2) ),
#else
      : paramsCorrectlyOrdered( _r1->getComposition() < _r2->getComposition() ),
#endif // defined(USE_ORIG_REACTANT_COMP_STRING)
        first(paramsCorrectlyOrdered ? _r1 : _r2),
        second(paramsCorrectlyOrdered ? _r2 : _r1),
        kConstant(0.0) {
    }

public:

	/**
	 * The rate constant
	 */
	double kConstant;


    /**
     * First cluster in reaction pair.
     * Reactant concentration guaranteed to be <= that of second cluster.
     */
    IReactant* first;

    /**
     * Second cluster in reaction pair.
     * Reactant concentration guaranteed to be >= that of first cluster.
     */
    IReactant* second;

    /**
     * Default constructor, deleted to ensure we are constructed with reactants.
     */
    Reaction() = delete;
};


/**
 * A Reaction with a descriptive key.
 * We separate out the descriptive key handling from the base Reaction
 * class for two reasons: we want a base Reaction class that is not
 * templatized, and we want our derived classes to be able to use
 * the key type that is most appropriate to their situation.
 */
template<typename K>
class KeyedReaction : public Reaction
{
public:
	/**
	 * Type of a canonical key describing this ProductionReaction that
	 * can be used to compare it to other ProductionReactions.
	 */
    using KeyType = K;
    

protected:
	/**
	 * A descriptive key in canonical form that can be used for
	 * fast compares against that of other Reactions of the same
     * derived type.  (No guarantees regarding comparisons between
     * keys for different types of reactions.)
	 */
	KeyType descKey;


    /**
     * Construct a KeyedReaction.
     * @param _r1 One of the reactants.
     * @param _r2 The other reactant.
     */
    KeyedReaction(IReactant* _r1, IReactant* _r2)
      : Reaction(_r1, _r2) {
    }

public:
    /**
     * Default and copy constructors, deleted to prevent use.
     */
    KeyedReaction() = delete;
    KeyedReaction(const KeyedReaction& other) = delete;

	/**
	 * Find the canonical key describing this ProductionReaction.
	 */
	const KeyType& descriptiveKey() const {
		return descKey;
	}
};

}

#endif /* XCORE_REACTION_H */
