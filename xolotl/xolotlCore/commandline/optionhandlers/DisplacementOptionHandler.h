#ifndef DISPLACEMENTOPTIONHANDLER_H
#define DISPLACEMENTOPTIONHANDLER_H

// Includes
#include <stdlib.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * DisplacementOptionHandler handles the displacement option.
 */
class DisplacementOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	DisplacementOptionHandler() :
    	OptionHandler("displacement",
    			"krFluence <value>                      "
    			"thresholdEnergy <value>                      "
    			"This option allows the user to change the Kr fluence "
    			"and displacement threshold energy.\n") {}

	/**
	 * The destructor
	 */
    ~DisplacementOptionHandler() {}

    /**
     * This method will set the IOptions krFluence and thresholdEnergy
     * to the value given as the argument.
     *
     * @param arg The value for the krypton fluence.
     * @param arg The value for the threshold energy.
     */
    bool handler(IOptions *opt, const std::string& arg) {

    	// Set the value for the fluence
    	double krFluence = strtod(arg.c_str(), NULL);

    	opt->setKrFluenceAmplitude(krFluence);

    	// Set the value for the threshold displacement energy
    	double thresholdEnergy = strtod(arg.c_str(), NULL);

    	opt->setDispEnergy(thresholdEnergy);

    	return true;
    }

};//end class DisplacementOptionHandler

} /* namespace xolotlCore */

#endif
