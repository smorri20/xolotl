#ifndef MEDIUMFACTORY_H
#define MEDIUMFACTORY_H

#include <memory>
#include <IMediumFactory.h>
#include <TokenizedLineReader.h>
#include <XGBAdvectionHandler.h>
#include <YGBAdvectionHandler.h>
#include <ZGBAdvectionHandler.h>

namespace xolotlFactory {

/**
 * Realizes the IMediumFactory interface. Handles the flux, advection, and diffusion
 * for a specific medium.
 */
class MediumFactory: public IMediumFactory {
protected:

	//! The flux handler
	std::shared_ptr<xolotlCore::IFluxHandler> theFluxHandler;

	//! The advection handler
	std::vector<std::shared_ptr<xolotlCore::IAdvectionHandler> > theAdvectionHandler;

	//! The diffusion handler
	std::shared_ptr<xolotlCore::IDiffusionHandler> theDiffusionHandler;

public:

	/**
	 * The constructor
	 */
	MediumFactory() {
	}

	/**
	 * The destructor
	 */
	~MediumFactory() {
	}

	/**
	 * Initialize the medium conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	void initializeMedium(xolotlCore::Options &options) {
		// Wrong if both he flux and time profile options are used
		if (options.useFluxAmplitude() && options.useFluxTimeProfile()) {
			// A constant flux value AND a time profile cannot both be given.
			throw std::string(
					"\nA constant flux value AND a time profile cannot both be given.");
		} else if (options.useFluxAmplitude()) {
			// Set the constant value of the flux
			theFluxHandler->setFluxAmplitude(options.getFluxAmplitude());
		} else if (options.useFluxTimeProfile()) {
			// Initialize the time profile
			theFluxHandler->initializeTimeProfile(options.getFluxProfileName());
		}

		// Get the number of dimensions
		int dim = options.getDimensionNumber();

		// Set-up the grain boundaries from the options
		std::string gbString = options.getGbString();
		// Build an input stream from the GB string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared < std::istringstream > (gbString);
		reader.setInputStream(argSS);
		// Break the string into tokens.
		auto tokens = reader.loadLine();
		// Loop on them
		for (int i = 0; i < tokens.size(); i++) {
			// Switch on the type of grain boundaries
			if (tokens[i] == "X") {
				auto GBAdvecHandler = std::make_shared<
						xolotlCore::XGBAdvectionHandler>();
				GBAdvecHandler->setPosition(
						strtod(tokens[i + 1].c_str(), NULL));
				theAdvectionHandler.push_back(GBAdvecHandler);
			} else if (tokens[i] == "Y") {
				if (dim < 2)
					// A Y grain boundary cannot be used in 1D.
					throw std::string(
							"\nA Y grain boundary CANNOT be used in 1D. Switch to 2D or 3D or remove it.");

				auto GBAdvecHandler = std::make_shared<
						xolotlCore::YGBAdvectionHandler>();
				GBAdvecHandler->setPosition(
						strtod(tokens[i + 1].c_str(), NULL));
				theAdvectionHandler.push_back(GBAdvecHandler);
			} else if (tokens[i] == "Z") {
				if (dim < 3)
					// A Z grain boundary cannot be used in 1D/2D.
					throw std::string(
							"\nA Z grain boundary CANNOT be used in 1D/2D. Switch to 3D or remove it.");

				auto GBAdvecHandler = std::make_shared<
						xolotlCore::ZGBAdvectionHandler>();
				GBAdvecHandler->setPosition(
						strtod(tokens[i + 1].c_str(), NULL));
				theAdvectionHandler.push_back(GBAdvecHandler);
			} else {
				// Wrong GB type
				throw std::string(
						"\nThe type of grain boundary is not known: \""
								+ tokens[i] + "\"");
			}

			i++;
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
	std::vector<std::shared_ptr<xolotlCore::IAdvectionHandler> > getAdvectionHandler() const {
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
};

} // end namespace xolotlFactory

#endif // MATERIALFACTORY_H
