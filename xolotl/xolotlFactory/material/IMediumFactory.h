#ifndef IMEDIUMFACTORY_H
#define IMEDIUMFACTORY_H

#include <memory>
#include <Options.h>
#include <IFluxHandler.h>
#include <IAdvectionHandler.h>
#include <IDiffusionHandler.h>

namespace xolotlFactory {

/**
 * Realizations of this interface are responsible for handling the flux, the advection,
 * and the diffusion. They are dependent on the medium under study.
 */
class IMediumFactory {
public:

	/**
	 * The destructor
	 */
	~IMediumFactory() {
	}

	/**
	 * Initialize the medium conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	virtual void initializeMedium(xolotlCore::Options &options) = 0;

	/**
	 * Return the flux handler.
	 *
	 * @return The flux handler.
	 */
	virtual std::shared_ptr<xolotlCore::IFluxHandler> getFluxHandler() const = 0;

	/**
	 * Return the advection handlers.
	 *
	 * @return The advection handlers.
	 */
	virtual std::vector<std::shared_ptr<xolotlCore::IAdvectionHandler> > getAdvectionHandler()
			const = 0;

	/**
	 * Return the diffusion handler.
	 *
	 * @return The diffusion handler.
	 */
	virtual std::shared_ptr<xolotlCore::IDiffusionHandler> getDiffusionHandler() const = 0;

	/**
	 * Function that creates and returns the wanted medium factory depending on the given
	 * orientation.
	 *
	 * @param orientation The orientation of the medium.
	 * @param dimension The number of dimensions of the problem.
	 * @return The medium factory.
	 */
	static std::shared_ptr<IMediumFactory> createMediumFactory(
			const std::string& orientation, int dimension);

};

} // end namespace xolotlFactory

#endif // IMEDIUMFACTORY_H
