#ifndef W100MEDIUMFACTORY_H
#define W100MEDIUMFACTORY_H

#include <memory>
#include <MediumFactory.h>
#include <W100FitFluxHandler.h>
#include <W100AdvectionHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MediumFactory for a (100) oriented tungsten medium.
 */
class W100MediumFactory : public MediumFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W100MediumFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W100MediumFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W100FitFluxHandler>();
		theAdvectionHandler.push_back(std::make_shared<xolotlCore::W100AdvectionHandler>());

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
				throw std::string("\nxolotlFactory: Bad dimension for the W100 medium factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W100MediumFactory() {}
};

} // end namespace xolotlFactory

#endif // W100MEDIUMFACTORY_H
