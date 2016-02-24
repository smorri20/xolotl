#ifndef WSRIMMATERIALHANDLERFACTORY_H
#define WSRIMMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <WSRIMFitFluxHandler.h>
#include <WSRIMAdvectionHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (SRIM) tungsten material.
 */
class WSRIMMaterialFactory : public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	WSRIMMaterialFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	WSRIMMaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::WSRIMFitFluxHandler>();
		theAdvectionHandler = std::make_shared<xolotlCore::WSRIMAdvectionHandler>();

		// Switch on the dimension for the diffusion handler
		switch (dim) {
			case 1:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion1DHandler>();
				break;
			case 2:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion2DHandler>();
				break;
			case 3:
				theDiffusionHandler = std::make_shared<xolotlCore::Diffusion3DHandler>();
				break;
			default:
				// The asked dimension is not good (e.g. -1, 0, 4)
				throw std::string("\nxolotlFactory: Bad dimension for the WSRIM material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~WSRIMMaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // WSRIMMATERIALHANDLERFACTORY_H
