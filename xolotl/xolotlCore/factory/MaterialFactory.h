#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include <IMaterialFactory.h>
#include <IFluxHandler.h>
#include <IAdvectionHandler.h>

namespace xolotlCore {

/**
 * Realizes the IMaterialFactory interface. Handles the flux and the advection
 * for a specific material.
 */
class MaterialFactory: public IMaterialFactory {
protected:

	//! The flux handler
	std::shared_ptr<IFluxHandler> theFluxHandler;

	//! The advection handler
	std::shared_ptr<IAdvectionHandler> theAdvectionHandler;

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
	void initializeMaterial(Options &options) {
		// If the Helium fluence option is present, set the value
		if (options.useMaxHeliumFluence()) {
			theFluxHandler->setMaxHeFluence(options.getMaxHeliumFluence());
		}

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

		return;
	}

	/**
	 * Return the flux handler.
	 *
	 *  @return The flux handler.
	 */
	std::shared_ptr<IFluxHandler> getFluxHandler() const {
		return theFluxHandler;
	}

	/**
	 * Return the advection handler.
	 *
	 *  @return The advection handler.
	 */
	std::shared_ptr<IAdvectionHandler> getAdvectionHandler() const {
		return theAdvectionHandler;
	}
};

} // end namespace xolotlCore

#endif // MATERIALHANDLERFACTORY_H
