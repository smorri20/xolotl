#include "DummyHandlerRegistry.h"

namespace xolotlPerf {

// Obtain a Timer by name.
std::shared_ptr<ITimer> DummyHandlerRegistry::getTimer(
		const std::string& name) {
	// TODO is there a need for us to retain access to this Timer?
	// TODO do we need to check whether client has already created
	// an object with this name and return that object?
	return std::make_shared < DummyTimer > (name);
}

// Obtain an EventCounter by name.
std::shared_ptr<IEventCounter> DummyHandlerRegistry::getEventCounter(
		const std::string& name) {
	// TODO is there a need for us to retain access to this Timer?
	// TODO do we need to check whether client has already created
	// an object with this name and return that object?
	return std::make_shared < DummyEventCounter > (name);
}

// Obtain a HardwareCounter object by name and by the
// counter data it collects.
std::shared_ptr<IHardwareCounter> DummyHandlerRegistry::getHardwareCounter(
		const std::string& name, const IHardwareCounter::SpecType& ctrSpec) {
	// TODO is there a need for us to retain access to this Timer?
	// TODO do we need to check whether client has already created
	// an object with this name and return that object?
	return std::make_shared < DummyHardwareCounter > (name, ctrSpec);
}

/**
 * Obtain a memory sampling region.
 */
std::shared_ptr<IMemSamplingRegion> DummyHandlerRegistry::getMemSamplingRegion(
        const std::string& name) {

    return std::make_shared<DummyMemSamplingRegion>(name);
}

} // namespace xolotlPerf

