#ifndef KRFLUENCEOPTIONHANDLER_H
#define KRFLUENCEOPTIONHANDLER_H

// Includes
#include <stdlib.h>
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * KrFluenceOptionHandler handles the displacement option.
 */
class KrFluenceOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	KrFluenceOptionHandler() :
    	OptionHandler("krFluence",
    			"krFluence <value>                      "
    			"This option allows the user to change the "
    			"Kr fluence .\n") {}

	/**
	 * The destructor
	 */
    ~KrFluenceOptionHandler() {}

    /**
     * This method will set the IOptions krFluence
     * to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The value for the krypton fluence.
     */
    bool handler(IOptions *opt, const std::string& arg) {

    	// Set the value for the fluence
    	double krFluence = strtod(arg.c_str(), NULL);

    	opt->setKrFluenceAmplitude(krFluence);

    	return true;
    }

};//end class KrFluenceOptionHandler

} /* namespace xolotlCore */

#endif
