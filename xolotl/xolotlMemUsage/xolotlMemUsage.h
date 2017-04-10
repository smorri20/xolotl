#ifndef XOLOTL_MEM_USAGE_H
#define XOLOTL_MEM_USAGE_H

#include <memory>
#include <sstream>
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
inline IHandlerRegistry::RegistryType toRegistryType(
		const std::string& arg) {

	IHandlerRegistry::RegistryType ret;

	if (arg == "dummy") {
		ret = IHandlerRegistry::dummy;
	} else if (arg == "std") {
		ret = IHandlerRegistry::std;
	} else {
		std::ostringstream estr;
		estr << "Invalid memory usage handler argument \"" << arg << "\" seen.";
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
                    std::chrono::duration<double>(1));

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
