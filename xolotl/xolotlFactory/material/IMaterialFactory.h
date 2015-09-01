#ifndef IMATERIALFACTORY_H
#define IMATERIALFACTORY_H

#include <memory>
#include <IMediumFactory.h>

namespace xolotlFactory {

/**
 * Realizations of this interface are responsible for creating different medium factories
 * that will represent the whole material under consideration.
 */
class IMaterialFactory {
public:

	/**
	 * The destructor
	 */
	~IMaterialFactory() {
	}

	/**
	 * Initialize the material conditions.
	 *
	 * @param options The Xolotl options
	 */
	virtual void initializeMaterial(xolotlCore::Options &options) = 0;

	/**
	 * Return the medium factory at a specific position.
	 *
	 * @param position The position on the grid
	 * @return The medium factory
	 */
	virtual std::shared_ptr<IMediumFactory> getMediumFactory(
			std::vector<double> position) const = 0;

	/**
	 * Return all the medium factories.
	 *
	 * @return The medium factories
	 */
	virtual std::vector<std::shared_ptr<IMediumFactory> > getMaterial() const = 0;

	/**
	 * Function that creates and returns the wanted material factory depending on the
	 * given options.
	 *
	 * @param options The Xolotl option
	 * @return The material factory
	 */
	static std::shared_ptr<IMaterialFactory> createMaterialFactory(xolotlCore::Options &options);

};

} // end namespace xolotlFactory

#endif // IMATERIALFACTORY_H
