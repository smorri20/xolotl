#include <iostream>
#include <sstream>
#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/dummy/DummyHandlerRegistry.h"
#include "xolotlMemUsage/standard/StdHandlerRegistry.h"


namespace xolotlMemUsage {

static std::shared_ptr<IHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
void initialize(IHandlerRegistry::RegistryType rtype,
                IHandlerRegistry::SamplingInterval samplingInterval) {

	switch (rtype) {
	case IHandlerRegistry::dummy:
		theHandlerRegistry = std::make_shared<DummyHandlerRegistry>(samplingInterval);
		break;

	case IHandlerRegistry::std:
        theHandlerRegistry = std::make_shared<StdHandlerRegistry>(samplingInterval);
		break;

	default:
		throw std::invalid_argument(
				"unrecognized memory usage handler registry type requested");
		break;
	}
}

// Provide access to our handler registry.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry(void) {
	if (!theHandlerRegistry) {
		throw std::runtime_error(
				"Request for xolotlMemUsage handler registry before xolotlMemUsage library has been initialized");
	}
	return theHandlerRegistry;
}

} // end namespace xolotlMemUsage

