#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include <IMaterialFactory.h>
#include <DummyAdvectionHandler.h>
#include <DummyDiffusionHandler.h>
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

	//! The diffusion handler
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

		// Check if all the physics processes should be used
		std::string process = options.getPhysicsProcess();

		// Don't do anything if all the processes are to be used
		if (process == "all") return;

		// Check if we just want the diffusion
		if (process == "diff") {
			// Turn off the advection
			theAdvectionHandler = std::make_shared<xolotlCore::DummyAdvectionHandler>();
			// Turn off the bubble bursting
			theBubbleBurstingHandler = std::make_shared<xolotlCore::DummyBubbleBurstingHandler>();
		}

		// Check if we just want the advection
		else if (process == "advec") {
			// Turn off the diffusion
			theDiffusionHandler = std::make_shared<xolotlCore::DummyDiffusionHandler>();
			// Turn off the bubble bursting
			theBubbleBurstingHandler = std::make_shared<xolotlCore::DummyBubbleBurstingHandler>();
		}

		// Check if we just want the bubble bursting
		else if (process == "burst") {
			// Turn off the diffusion
			theDiffusionHandler = std::make_shared<xolotlCore::DummyDiffusionHandler>();
			// Turn off the advection
			theAdvectionHandler = std::make_shared<xolotlCore::DummyAdvectionHandler>();
		}

		// Check if we don't want any of them
		else if (process == "none") {
			// Turn off the diffusion
			theDiffusionHandler = std::make_shared<xolotlCore::DummyDiffusionHandler>();
			// Turn off the advection
			theAdvectionHandler = std::make_shared<xolotlCore::DummyAdvectionHandler>();
			// Turn off the bubble bursting
			theBubbleBurstingHandler = std::make_shared<xolotlCore::DummyBubbleBurstingHandler>();
		}

		// Wrong option name
		else {
			throw std::string(
							"\nThe physics process name is not known: \"" + process
							+ "\", look at Xolotl's help to know which names are available.");
		}

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
	 * Return the bubble bursting handler.
	 *
	 * @return The bubble bursting handler.
	 */
	std::shared_ptr<xolotlCore::IBubbleBurstingHandler> getBubbleBurstingHandler() const {
		return theBubbleBurstingHandler;
	}
};

} // end namespace xolotlFactory

#endif // MATERIALHANDLERFACTORY_H
