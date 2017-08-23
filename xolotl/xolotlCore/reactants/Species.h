#ifndef XOLOTL_CORE_SPECIES_H_
#define XOLOTL_CORE_SPECIES_H_

#include <string>

namespace xolotlCore {


enum class Species {
    Invalid = -1,
    H,
    He,
    V,
    I,
    HeV,
    HeI,
    Xe,
    XeV,
    XeI,
    NESuper,
    PSISuper
};

std::string toString(const Species& s);

} // namespace xolotlCore

#endif /* XOLOTL_CORE_SPECIES_H_ */
