#ifndef W111MATERIALHANDLERFACTORY_H
#define W111MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W111FitFluxHandler.h>
#include <W111AdvectionHandler.h>
#include <Diffusion1DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (11) oriented tungsten material.
 */
class W111MaterialFactory : public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 */
	W111MaterialFactory() {
		theFluxHandler = std::make_shared<xolotlCore::W111FitFluxHandler>();
		theAdvectionHandler = std::make_shared<xolotlCore::W111AdvectionHandler>();
		theDiffusionHandler = std::make_shared<xolotlCore::Diffusion1DHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W111MaterialFactory() {}
};

} // end namespace xolotlFactory

#endif // W111MATERIALHANDLERFACTORY_H
