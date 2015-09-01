#include <IMaterialFactory.h>
#include <MaterialFactory.h>

namespace xolotlFactory {

static std::shared_ptr<IMaterialFactory> theMaterialFactory;

std::shared_ptr<IMaterialFactory> IMaterialFactory::createMaterialFactory(xolotlCore::Options &options) {
	theMaterialFactory = std::make_shared<MaterialFactory>(options);

	return theMaterialFactory;
}

};  // end namespace xolotlFactory
