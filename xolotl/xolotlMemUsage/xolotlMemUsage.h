#ifndef XOLOTL_MEM_USAGE_H
#define XOLOTL_MEM_USAGE_H

#include <memory>
#include <sstream>
#include <iterator>
#include "xolotlMemUsage/IHandlerRegistry.h"
#include "xolotlMemUsage/RuntimeError.h"

namespace xolotlMemUsage {

/**
 * Detect the type of handlers registry to create based
 * on a string argument (e.g., taken from the command line).
 * Throws an std::invalid_argument exception if arg does not
 * specify a registry type we know about.
 *
 * NOTE: We define the function in the header to avoid adding a dependency
 * on the library for the command line library which uses
 * this function to parse its command line.
 *
 * @param arg Type of handler registry to create.
 * @return Newly-created handler registry.
 */
inline IHandlerRegistry::RegistryType toRegistryType(const std::vector<std::string>& tokens) {

	IHandlerRegistry::RegistryType ret;

    try {
        if(tokens[0] == "dummy") {
            ret = IHandlerRegistry::dummy;
        } else if(tokens[0] == "summary") {
            if(tokens[1] == "proc") {
                ret = IHandlerRegistry::summaryProc;
            } else if(tokens[1] == "node") {
                ret = IHandlerRegistry::summaryNode;
            } else {
                throw std::invalid_argument("unrecognized registry variant");
            }
        } else if(tokens[0] == "profile") {
            if(tokens[1] == "proc") {
                ret = IHandlerRegistry::profileProc;
            } else if(tokens[1] == "node") {
                ret = IHandlerRegistry::profileNode;
            } else {
                throw std::invalid_argument("unrecognized registry variant");
            }
        } else  {
            throw std::invalid_argument("unrecognized primary registry type");
        }
    }
    catch(std::exception& e) {
        std::ostringstream estr;
        estr << "Invalid memory usage handler arguments: ";
        std::copy(tokens.begin(), tokens.end(),
                std::ostream_iterator<std::string>(estr, " "));
        estr << '\n' << e.what() << '\n';
		throw std::invalid_argument(estr.str());
    }
	return ret;
}

/**
 * Initialize the memory usage library for using the desired type of handlers.
 * Throws a std::invalid_argument if caller requests a registry type that
 * we do not support.
 *
 * @param rtype Type of handlerRegistry to create.
 * @param samplingInterval Duration between samples for asynchronously-sampled data (default 1s).
 */
void initialize(IHandlerRegistry::RegistryType rtype,
                IHandlerRegistry::SamplingInterval samplingInterval = 
                    std::chrono::duration<uint64_t>(1),
                std::string profileFilename = "");

/**
 * Access the handler registry.
 * Throws a std::runtime_error if called before the xolotlMemUsage classes
 * have been initialized.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<IHandlerRegistry> getHandlerRegistry(void);

} // end namespace xolotlMemUsage

#endif // XOLOTL_MEM_USAGE_H
