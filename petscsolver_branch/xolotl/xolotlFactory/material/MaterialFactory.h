#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include <IMaterialFactory.h>
#include <DummyAdvectionHandler.h>
#include <DummyDiffusionHandler.h>
#include <DummyTrapMutationHandler.h>
#include <DummyBubbleBurstingHandler.h>

namespace xolotlFactory {

/**
 * Realizes the IMaterialFactory interface. Handles the flux and the advection
 * for a specific material.
 */
class MaterialFactory: public IMaterialFactory {
protected:

	//! The flux handler
	std::shared_ptr<xolotlCore::IFluxHandler> theFluxHandler;

	//! The advection handler
	std::shared_ptr<xolotlCore::IAdvectionHandler> theAdvectionHandler;

	//! The diffusion handler
	std::shared_ptr<xolotlCore::IDiffusionHandler> theDiffusionHandler;

	//! The modified trap-mutation handler
	std::shared_ptr<xolotlCore::ITrapMutationHandler> theTrapMutationHandler;

	//! The bubble bursting handler
	std::shared_ptr<xolotlCore::IBubbleBurstingHandler> theBubbleBurstingHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	MaterialFactory() {
	}

	/**
	 * The destructor
	 */
	~MaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	void initializeMaterial(xolotlCore::Options &options) {
		// Wrong if both he flux and time profile options are used
		if (options.useHeliumFlux() && options.useFluxTimeProfile()) {
			// A constant flux value AND a time profile cannot both be given.
			throw std::string(
					"\nA constant flux value AND a time profile cannot both be given.");
		}
		else if (options.useHeliumFlux()) {
			// Set the constant value of the flux
			theFluxHandler->setHeFlux(options.getHeliumFlux());
		}
		else if (options.useFluxTimeProfile()) {
			// Initialize the time profile
			theFluxHandler->initializeTimeProfile(options.getFluxProfileName());
		}

		// Get the process map
		auto map = options.getProcesses();
		// Set dummy handlers when needed
		if (!map["diff"])
			theDiffusionHandler = std::make_shared<xolotlCore::DummyDiffusionHandler>();
		if (!map["advec"])
			theAdvectionHandler = std::make_shared<xolotlCore::DummyAdvectionHandler>();
		if (!map["modifiedTM"])
			theTrapMutationHandler = std::make_shared<xolotlCore::DummyTrapMutationHandler>();
		if (!map["bursting"])
			theBubbleBurstingHandler = std::make_shared<xolotlCore::DummyBubbleBurstingHandler>();

		return;
	}

	/**
	 * Return the flux handler.
	 *
	 *  @return The flux handler.
	 */
	std::shared_ptr<xolotlCore::IFluxHandler> getFluxHandler() const {
		return theFluxHandler;
	}

	/**
	 * Return the advection handler.
	 *
	 *  @return The advection handler.
	 */
	std::shared_ptr<xolotlCore::IAdvectionHandler> getAdvectionHandler() const {
		return theAdvectionHandler;
	}

	/**
	 * Return the diffusion handler.
	 *
	 *  @return The diffusion handler.
	 */
	std::shared_ptr<xolotlCore::IDiffusionHandler> getDiffusionHandler() const {
		return theDiffusionHandler;
	}

	/**
	 * Return the modified trap-mutation handler.
	 *
	 *  @return The trap-mutation handler.
	 */
	std::shared_ptr<xolotlCore::ITrapMutationHandler> getTrapMutationHandler() const {
		return theTrapMutationHandler;
	}

	/**
	 * Return the bubble bursting handler.
	 *
	 *  @return The bubble bursting handler.
	 */
	std::shared_ptr<xolotlCore::IBubbleBurstingHandler> getBubbleBurstingHandler() const {
		return theBubbleBurstingHandler;
	}
};

} // end namespace xolotlFactory

#endif // MATERIALHANDLERFACTORY_H
