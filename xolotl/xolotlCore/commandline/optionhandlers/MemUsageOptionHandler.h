#ifndef MEMUSAGEOPTIONHANDLER_H
#define MEMUSAGEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"
#include "xolotlMemUsage/IHandlerRegistry.h"

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
				"memUsageHandler dummy | std [sampling_interval_ms] | (profileproc|profilenode) output_filename [sampling_interval_ms]\n"
				"Which memory usage handlers to use (default = std) and sampling interval in ms (default=1000, ignored with dummy handler)\n") {}

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
            // Tokenize the input line.
		    xolotlCore::TokenizedLineReader<std::string> tokenizer;
            tokenizer.setInputStream(std::make_shared<std::istringstream>(arg));
            auto tokens = tokenizer.loadLine();

            if(tokens.size() < 1)
            {
                throw std::invalid_argument("invalid memUsageHandler options given.");
            }

            // Determine type of handler to use.
            xolotlMemUsage::IHandlerRegistry::RegistryType rtype = 
                xolotlMemUsage::toRegistryType(tokens[0]);
            opt->setMemUsageHandlerType(rtype);

            // Handle additional arguments if necessary.
            if(rtype != xolotlMemUsage::IHandlerRegistry::dummy) {
                
                auto currIdx = 1;

                if((rtype == xolotlMemUsage::IHandlerRegistry::profileproc) or
                    (rtype == xolotlMemUsage::IHandlerRegistry::profilenode)) {

                    if(tokens.size() < 2) {
                        // We were not given a filename to use.
                        throw std::invalid_argument("Memory usage profiling requested, but no output filename given for memory usage profiles");
                    }

                    opt->setMemUsageProfileFilename(tokens[currIdx]);
                    ++currIdx;
                }
                
                // See if we were given a sampling interval.
                if(currIdx < tokens.size()) {
                    std::istringstream sistr(tokens[currIdx]);
                    uint64_t samplingIntervalMillis = 0;
                    sistr >> samplingIntervalMillis;
                    if(not sistr.eof())
                    {
                        throw std::invalid_argument("unable to convert mem usage sampling interval argument to ms");
                    }
                    if(samplingIntervalMillis == 0)
                    {
                        throw std::invalid_argument("sampling interval must be > 0 ms");
                    }
                    opt->setMemUsageSamplingInterval(std::chrono::duration<uint64_t, std::milli>(samplingIntervalMillis));
                }
            }
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
