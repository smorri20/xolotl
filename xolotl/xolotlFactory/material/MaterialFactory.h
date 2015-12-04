#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include <IMaterialFactory.h>
#include <Td40FitDisplacementHandler.h>
#include <W100FitFluxHandler.h>

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

	//! The displacement handler
	std::shared_ptr<xolotlCore::IDisplacementHandler> theDisplacementHandler;


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
		if (options.useFluxAmplitude() && options.useFluxTimeProfile()) {
			// A constant flux value AND a time profile cannot both be given.
			throw std::string(
					"\nA constant flux value AND a time profile cannot both be given.");
		}
		else if (options.useFluxAmplitude()) {
			// Set the constant value of the flux
			theFluxHandler->setFluxAmplitude(options.getFluxAmplitude());
		}
		else if (options.useFluxTimeProfile()) {
			// Initialize the time profile
			theFluxHandler->initializeTimeProfile(options.getFluxProfileName());
		}

//		theDisplacementHandler->setKrFluenceAmplitude(options.getKrFluenceAmplitude());
//		theDisplacementHandler->setDispEnergy(options.getDispEnergy());

//		theFluxHandler = std::make_shared<xolotlCore::W100FitFluxHandler>();
    	theDisplacementHandler = std::make_shared<xolotlCore::Td40FitDisplacementHandler>();

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
	 * Return the displacement handler.
	 *
	 *  @return The displacement handler.
	 */
	std::shared_ptr<xolotlCore::IDisplacementHandler> getDisplacementHandler() const {
		return theDisplacementHandler;
	}

};

} // end namespace xolotlFactory

#endif // MATERIALHANDLERFACTORY_H
