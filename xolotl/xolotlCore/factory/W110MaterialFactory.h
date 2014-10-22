#ifndef W110MATERIALHANDLERFACTORY_H
#define W110MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W110FitFluxHandler.h>
#include <W110AdvectionHandler.h>

namespace xolotlCore {

/**
 * Subclass of MaterialFactory for a (110) oriented tungsten material.
 */
class W110MaterialFactory : public MaterialFactory {

public:

	/**
	 * The constructor creates the handlers.
	 */
	W110MaterialFactory() {
		theFluxHandler = std::make_shared<W110FitFluxHandler>();
		theAdvectionHandler = std::make_shared<W110AdvectionHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W110MaterialFactory() {}
};

} // end namespace xolotlCore

#endif // W110MATERIALHANDLERFACTORY_H
