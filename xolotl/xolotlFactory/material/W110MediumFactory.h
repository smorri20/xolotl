#ifndef W110MEDIUMFACTORY_H
#define W110MEDIUMFACTORY_H

#include <memory>
#include <MediumFactory.h>
#include <W110FitFluxHandler.h>
#include <W110AdvectionHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MediumFactory for a (110) oriented tungsten medium.
 */
class W110MediumFactory : public MediumFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W110MediumFactory() {}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W110MediumFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W110FitFluxHandler>();
		theAdvectionHandler.push_back(std::make_shared<xolotlCore::W110AdvectionHandler>());

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
				throw std::string("\nxolotlFactory: Bad dimension for the W110 medium factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W110MediumFactory() {}
};

} // end namespace xolotlFactory

#endif // W110MEDIUMFACTORY_H
