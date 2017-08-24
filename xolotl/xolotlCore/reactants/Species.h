#ifndef XOLOTL_CORE_SPECIES_H_
#define XOLOTL_CORE_SPECIES_H_

#include <string>

namespace xolotlCore {

/**
 * Reactant species that Xolotl supports.
 * Note that these are individual species, *not* cluster types (so we don't 
 * include things like HeV here).
 * Note also that the enum values for valid species should be
 * consecutive integers starting with 0, so that we can use them
 * as indices into an IReactant::Composition array.
 */
/*
 *
 * TODO if we change Xolotl configuration so it builds only one type of
 * reaction network, this should go into the *clusters directory so that
 * it is specific to the type of reaction network used.
 */
enum class Species {
    Invalid = -1,
    V = 0,
    I,
    He,
    Xe,
    last = Xe
};

/**
 * Obtain a human-readable string name for the given Species.
 *
 * @param s The Species of interest.
 * @return A string representation of the given Species s.
 */
std::string toString(Species s);


/**
 * Obtain the integer value of a Species.
 * We define this so that conversions are explicit and their intended use
 * is made clear in the code that uses the conversion.
 *
 * @param s The Species of interest.
 * @return The index within an IReactant::Composition of the concentration
 * of Species s in that composition.
 */
inline
int toCompIdx(Species s)
{
    return static_cast<int>(s);
}

} // namespace xolotlCore

#endif /* XOLOTL_CORE_SPECIES_H_ */
