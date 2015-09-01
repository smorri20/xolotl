#include <IMediumFactory.h>
#include <W100MediumFactory.h>
#include <W110MediumFactory.h>
#include <W111MediumFactory.h>

namespace xolotlFactory {

static std::shared_ptr<IMediumFactory> theMediumFactory;

std::shared_ptr<IMediumFactory> IMediumFactory::createMediumFactory(
		const std::string& orientation, int dimension) {
	// W100 case
	if (orientation == "W100")
		theMediumFactory = std::make_shared < W100MediumFactory > (dimension);
	// W110 case
	else if (orientation == "W110")
		theMediumFactory = std::make_shared < W110MediumFactory > (dimension);
	// W111 case
	else if (orientation == "W111")
		theMediumFactory = std::make_shared < W111MediumFactory > (dimension);
	// The type is not supported
	else {
		throw std::string(
				"\nThe material type is not known: \"" + orientation + "\"");
	}
	return theMediumFactory;
}

}
;
// end namespace xolotlFactory
