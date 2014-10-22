#ifndef W111MATERIALHANDLERFACTORY_H
#define W111MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W111FitFluxHandler.h>
#include <W111AdvectionHandler.h>

namespace xolotlCore {

/**
 * Subclass of MaterialFactory for a (100) oriented tungsten material.
 */
class W111MaterialFactory : public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 */
	W111MaterialFactory() {
		theFluxHandler = std::make_shared<W111FitFluxHandler>();
		theAdvectionHandler = std::make_shared<W111AdvectionHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W111MaterialFactory() {}
};

} // end namespace xolotlCore

#endif // W111MATERIALHANDLERFACTORY_H
