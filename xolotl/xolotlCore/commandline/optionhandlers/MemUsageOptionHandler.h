#ifndef MEMUSAGEOPTIONHANDLER_H
#define MEMUSAGEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * MemUsageOptionHandler handles the choice of handlers for the 
 * memory usage infrastructure.
 */
class MemUsageOptionHandler: public OptionHandler {
public:

	/**
	 * Construct a MemUsageOptionHandler.
	 */
	MemUsageOptionHandler() :
		OptionHandler("memUsageHandler",
				"memUsageHandler {std,dummy}   "
				"Which memory usage handlers to use. (default = std)\n") {}

	/**
	 * Destroy the MemUsageOptionHandler.
	 */
	~MemUsageOptionHandler() {
	}

	/**
	 * This method will set the IOptions memUsageStandardHandlersFlag
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The argument for the flag.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
        
        bool ret = true;

        try
        {
		    // Determine the type of handlers we are being asked to use
            xolotlMemUsage::IHandlerRegistry::RegistryType rtype = xolotlMemUsage::toRegistryType(arg);
            opt->setMemUsageHandlerType( rtype );
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << e.what() << std::endl;
            opt->showHelp(std::cerr);
            opt->setShouldRunFlag(false);
            opt->setExitCode(EXIT_FAILURE);
            ret = false;
        }

		return ret;
	}

};
//end class MemUsageOptionHandler

} /* namespace xolotlCore */

#endif
