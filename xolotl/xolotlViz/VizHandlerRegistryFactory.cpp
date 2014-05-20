#include "XolotlConfigViz.h"
#include "VizHandlerRegistryFactory.h"
#include <DummyHandlerRegistry.h>
#include <iostream>

#if defined(HAVE_VIZLIB_STD)
#include <StandardHandlerRegistry.h>
#endif // defined(HAVE_VIZLIB_STD)


namespace xolotlViz
{

static std::shared_ptr<IVizHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
void initialize(bool useStdRegistry) {

	if (useStdRegistry) {
#if defined(HAVE_VIZLIB_STD)
		// we are to use a standard handler registry
		theHandlerRegistry = std::make_shared<StandardHandlerRegistry>();
#else
		// the dependencies are missing, abort
		throw("xolotlViz::initialize: unable to build requested standard handler registry because the dependencies are missing.");
#endif // defined(HAVE_VIZLIB_STD)
	}
    else {
    	theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
    }

    return;
}

// Provide access to our handler registry.
std::shared_ptr<IVizHandlerRegistry> getVizHandlerRegistry( void )
{
    if( !theHandlerRegistry )
    {
    	// It is not initialized yet so throw an error
        throw("xolotlViz::getVizHandlerRegistry: VizHandlerRegistry has not been correctly initialized.");
    }

    return theHandlerRegistry;
}

} // end namespace xolotlViz

