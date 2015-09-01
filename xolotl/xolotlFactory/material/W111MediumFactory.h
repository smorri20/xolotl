#ifndef W111MEDIUMFACTORY_H
#define W111MEDIUMFACTORY_H

#include <memory>
#include <MediumFactory.h>
#include <W111FitFluxHandler.h>
#include <W111AdvectionHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MediumFactory for a (111) oriented tungsten medium.
 */
class W111MediumFactory : public MediumFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W111MediumFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W111MediumFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W111FitFluxHandler>();
		theAdvectionHandler.push_back(std::make_shared<xolotlCore::W111AdvectionHandler>());

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
				throw std::string("\nxolotlFactory: Bad dimension for the W111 medium factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W111MediumFactory() {}
};

} // end namespace xolotlFactory

#endif // W111MEDIUMFACTORY_H
