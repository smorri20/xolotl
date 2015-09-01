#ifndef MATERIALFACTORY_H
#define MATERIALFACTORY_H

#include <memory>
#include <IMaterialFactory.h>
#include <W100MediumFactory.h>
#include <W110MediumFactory.h>
#include <W111MediumFactory.h>

namespace xolotlFactory {

/**
 * Realizes the IMaterialFactory interface.
 */
class MaterialFactory: public IMaterialFactory {
private:

	// The vector of medium factories
	std::vector<std::shared_ptr<xolotlFactory::IMediumFactory> > theMediumFactories;

	// The vector of Y separations
	std::vector<double> theYSeparations;

	/**
	 * The constructor creates the handlers.
	 */
	MaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 */
	MaterialFactory(xolotlCore::Options &options) {
		// Set-up the material from the options
		std::string materialString = options.getMaterial();
		// Build an input stream from the material string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared < std::istringstream > (materialString);
		reader.setInputStream(argSS);
		// Break the string into tokens.
		auto tokens = reader.loadLine();
		// Get the total number of tokens
		int tokenSize = tokens.size();

		// The option is wrong if the number of token is even
		if (tokenSize % 2 == 0) {
			throw std::string(
					"\nThe format for the material is wrong, there is an even number of tokens!");
		}

		// There should be at list one material listed
		theMediumFactories.push_back(
				IMediumFactory::createMediumFactory(tokens[0],
						options.getDimensionNumber()));
		// Loop on the rest of the token
		for (int i = 1; i < tokens.size(); i++) {
			// Get where is the separation
			theYSeparations.push_back(strtod(tokens[i].c_str(), NULL));
			// Create the next medium
			theMediumFactories.push_back(
					IMediumFactory::createMediumFactory(tokens[i + 1],
							options.getDimensionNumber()));
			// Increment i
			i++;
		}

		return;
	}

	/**
	 * The destructor
	 */
	~MaterialFactory() {
	}

	/**
	 * Initialize the material conditions.
	 *
	 * @param options The Xolotl options
	 */
	void initializeMaterial(xolotlCore::Options &options) {
		// Loop on all the medium factories
		for (int i = 0; i < theMediumFactories.size(); i++) {
			theMediumFactories[i]->initializeMedium(options);
		}

		return;
	}

	/**
	 * Return the medium factory at a specific position.
	 *
	 * @return The medium factory
	 */
	std::shared_ptr<xolotlFactory::IMediumFactory> getMediumFactory(
			std::vector<double> position) const {
		// Return the first one if there is only one factory
		if (theMediumFactories.size() == 1) return theMediumFactories[0];

		// Return the first one if the Y position is lower than the first
		// stored separation
		if (position[1] < theYSeparations[0]) return theMediumFactories[0];

		// Find where the Y position is in the Y separations
		int indice = 0;
		for (int i = 0; i < theYSeparations.size(); i++) {
			if (position[1] < theYSeparations[i]) break;
			indice++;
		}

		// Return the corresponding medium factory
		return theMediumFactories[indice];
	}

	/**
	 * Return all the medium factories.
	 *
	 * @return The medium factories.
	 */
	std::vector<std::shared_ptr<IMediumFactory> > getMaterial() const {
		return theMediumFactories;
	}

}
;

} // end namespace xolotlFactory

#endif // MATERIALFACTORY_H
