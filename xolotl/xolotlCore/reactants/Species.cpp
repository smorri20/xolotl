#include <string>
#include <map>
#include "Species.h"

namespace xolotlCore {

std::string toString(const Species& s) {

    static std::map<Species, std::string> smap {
        { Species::Invalid, "Invalid_species" },
        { Species::H, "H" },
        { Species::He, "He" },
        { Species::V, "V" },
        { Species::I, "I" },
        { Species::HeV, "HeV" },
        { Species::HeI, "HeI" },
        { Species::Xe, "Xe" },
        { Species::XeV, "XeV" },
        { Species::XeI, "XeI" },
        { Species::NESuper, "NESuper" },
        { Species::PSISuper, "PSISuper" }
    };

    auto iter = smap.find(s);
    return (iter != smap.end()) ? iter->second : "[unrecognized species]";
}


} // namespace xolotlCore

