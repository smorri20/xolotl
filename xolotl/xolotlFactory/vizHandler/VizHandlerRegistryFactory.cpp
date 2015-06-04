#include <XolotlConfigViz.h>
#include <VizHandlerRegistryFactory.h>
#include <DummyHandlerRegistry.h>
#include <iostream>
#include <mpi.h>

#if defined(HAVE_VIZLIB_STD)
#include <StandardHandlerRegistry.h>
#endif // defined(HAVE_VIZLIB_STD)

namespace xolotlFactory {

static std::shared_ptr<xolotlViz::IVizHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
bool initializeVizHandler(xolotlViz::IVizHandlerRegistry::RegistryType rtype) {
	bool ret = true;

    switch ( rtype )
    {
      case xolotlViz::IVizHandlerRegistry::dummy:
		// we are to use a dummy handler registry
		theHandlerRegistry = std::make_shared<xolotlViz::DummyHandlerRegistry>();
        break;

#if defined(HAVE_VIZLIB_STD)

      case xolotlViz::IVizHandlerRegistry::std:
#if defined (HAVE_OSMESA)
        theHandlerRegistry = std::make_shared<xolotlViz::StandardHandlerRegistry>(
            xolotlViz::StandardHandlerRegistry::png);
#else
        theHandlerRegistry = std::make_shared<xolotlViz::StandardHandlerRegistry>(
            xolotlViz::StandardHandlerRegistry::eps);
#endif // defined(HAVE_OSMESA)
        break;


      case xolotlViz::IVizHandlerRegistry::png:
#if defined (HAVE_OSMESA)
        theHandlerRegistry = std::make_shared<xolotlViz::StandardHandlerRegistry>(
            xolotlViz::StandardHandlerRegistry::png);
#else
			throw std::string(
					"\nxolotlFactory::initialize: unable to build requested visualization "
							"png handler registry due to missing dependencies");
#endif // defined(HAVE_OSMESA)
        break;


      case xolotlViz::IVizHandlerRegistry::eps:
        theHandlerRegistry = std::make_shared<xolotlViz::StandardHandlerRegistry>(
            xolotlViz::StandardHandlerRegistry::eps);
        break;

#else // defined(HAVE_VIZLIB_STD)

      default:
      {
		// Get the current process ID
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		// Only print the error message once when running in parallel
		if (procId == 0) {
			// it is not possible to use the standard registry
			throw std::string(
					"\nxolotlFactory::initialize: unable to build requested visualization "
							"handler registry due to missing dependencies");
		}
      }
#endif // defined(HAVE_VIZLIB_STD)

    }

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<xolotlViz::IVizHandlerRegistry> getVizHandlerRegistry(void) {
	if (!theHandlerRegistry) {
		// Throw an error since we have not yet been initialized
		throw std::string(
				"\nxolotlFactory handler registry requested, but library has not been initialized.");
	}

	return theHandlerRegistry;
}

} // end namespace xolotlFactory
