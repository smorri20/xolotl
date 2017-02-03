#ifndef REACTANTUTILS_H
#define REACTANTUTILS_H

#include "IReactant.h"

namespace xolotlCore {

class IReactant;

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k is stored here and other information is implemented by
 * the daughter classes. k is computed when setTemperature() is called.
 */
class Reaction {
public:

	/**
	 * The rate constant
	 */
	double kConstant;

	//! The constructor
	Reaction() :
			kConstant(0.0) {
	}
};

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k+ is stored along the clusters taking part in the
 * production for faster computation because they only change
 * when the temperature change.
 * We don't keep the product cluster because it is not needed for k calculations.
 *
 * k is computed when setTemperature() is called.
 */
class ProductionReaction: public Reaction {
public:

	/**
	 * The first cluster in the pair
	 */
	IReactant * first;

	/**
	 * The second cluster in the pair
	 */
	IReactant * second;

	//! The constructor
	ProductionReaction(IReactant * firstPtr, IReactant * secondPtr) :
			first(firstPtr), second(secondPtr) {
	}
};

/**
 * This is a public class that is used to store a reaction.
 *
 * The constant k- is stored along the clusters taking part in the
 * production for faster computation because they only change
 * when the temperature change. k is computed when setTemperature() is called.
 */
class DissociationReaction: public Reaction {
public:

	/**
	 * The dissociating cluster
	 */
	IReactant * dissociating;

	/**
	 * The first cluster in the pair
	 */
	IReactant * first;

	/**
	 * The second cluster in the pair
	 */
	IReactant * second;

	/**
	 * The revert reaction corresponding to this dissociation
	 */
	ProductionReaction * reverseReaction;

	//! The constructor
	DissociationReaction(IReactant * dissociatingPtr, IReactant * firstPtr,
			IReactant * secondPtr) :
			dissociating(dissociatingPtr), first(firstPtr), second(secondPtr), reverseReaction(
					nullptr) {
	}
};

}

#endif /* REACTANTUTILS_H */
