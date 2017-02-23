#ifndef DISPLACEMENTOPTIONHANDLER_H
#define DISPLACEMENTOPTIONHANDLER_H

// Includes
#include <stdlib.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * DisplacementOptionHandler handles the displacement threshold energy option.
 */
class DisplacementOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	DisplacementOptionHandler() :
			OptionHandler("thresholdEnergy",
					"thresholdEnergy <value>                      "
							"This option allows the user to change the "
							"displacement threshold energy.\n") {
	}

	/**
	 * The destructor
	 */
	~DisplacementOptionHandler() {
	}

	/**
	 * This method will set the IOptions thresholdEnergy
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The value for the displacement threshold energy.
	 */
	bool handler(IOptions *opt, const std::string& arg) {

		// Set the value for the threshold displacement energy
		int thresholdEnergy = strtod(arg.c_str(), NULL);

		opt->setDispEnergy(thresholdEnergy);

		return true;
	}

};
//end class DisplacementOptionHandler

} /* namespace xolotlCore */

#endif
