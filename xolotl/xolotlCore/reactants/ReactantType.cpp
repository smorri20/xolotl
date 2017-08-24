#include <string>
#include <unordered_map>
#include "ReactantType.h"

namespace xolotlCore {

std::string toString(ReactantType rtype) {

    static std::unordered_map<ReactantType, std::string> smap {
        { ReactantType::Invalid, "Invalid_reactant_type" },
        { ReactantType::V, "V" },
        { ReactantType::I, "I" },
        { ReactantType::He, "He" },
        { ReactantType::HeV, "HeV" },
        { ReactantType::HeI, "HeI" },
        { ReactantType::PSISuper, "PSISuper" },
        { ReactantType::Xe, "Xe" },
        { ReactantType::XeV, "XeV" },
        { ReactantType::XeI, "XeI" },
        { ReactantType::NESuper, "NESuper" }
    };

    auto iter = smap.find(rtype);
    return (iter != smap.end()) ? iter->second : "[unrecognized reactant type]";
}


Species toSpecies(ReactantType rtype) {
    // Map of all single species reactant types, plus Invalid.
    static std::unordered_map<ReactantType, Species> smap {
        { ReactantType::Invalid, Species::Invalid },
        { ReactantType::V, Species::V },
        { ReactantType::I, Species::I },
        { ReactantType::He, Species::He },
        { ReactantType::Xe, Species::Xe }
    };

    auto iter = smap.find(rtype);
    return (iter != smap.end()) ? iter->second : Species::Invalid;
}


} // namespace xolotlCore

