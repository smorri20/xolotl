#include "HandlerRegistryFactory.h"
#include "dummy/DummyHandlerRegistry.h"

#if READY
#if have GPTL
// TODO this is the wrong test - we need to test whether 
// to build the standard classes, and then include the standard classes
// regardless of what they are built on top of.
#include "standard/StandardHandlerRegistry.h"
#endif // have GPTL

#else
#include "standard/StandardHandlerRegistry.h"
#endif // READY

namespace xolotlPerf
{

static std::shared_ptr<IHandlerRegistry> theHandlerRegistry;

// Parse the command line to figure out which type of handlerRegistry
// to create.  Then create the appropriate type of handlerRegistry
//
// we need something that will parse the command line and know
// what to build.  we don't want to require the user to have 
// GPTL and PAPI, and we want to have the "dummy" option to avoid
// taking times.
// So the dummy classes need to be included in every executable,
// and the GPTL-based ones need to be included.
void
initialize( int argc, char* argv[] )
{
#if READY
    // Parse the command line to determine which type of handler registry to make

    // Create the desired type of handler registry.
#else
    theHandlerRegistry = std::make_shared<StandardHandlerRegistry>();
#endif // READY
}



// Provide access to our handler registry.
std::shared_ptr<IHandlerRegistry>
getHandlerRegistry( void )
{
    return theHandlerRegistry;
}


} // end namedpsace xolotlPerf

