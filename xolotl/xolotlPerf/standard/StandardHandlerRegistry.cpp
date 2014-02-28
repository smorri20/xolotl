#include "StandardHandlerRegistry.h"
#include "gptl.h"


namespace xolotlPerf
{

#if RREADY
// currently commented out so that other recent modifications can be
// committed to the Subversion repository.

std::shared_ptr<StandardHandlerRegistry> theRegistry;


std::shared_ptr<StandardHandlerRegistry>
StandardHandlerRegistry::getRegistry( void )
{
    if( !theRegistry )
    {
        theRegistry = std::make_shared<StandardHandlerRegistry>;
    }
    return theRegistry;
}


StandardHandlerRegistry::StandardHandlerRegistry( void )
  : Identifiable( "theStandardHandlerRegistry" )
{
    // We use GPTL for data collection, so we must make sure 
    // that library has been initialized.
    //
    // Note: We assume that no other part of the code is using GPTL.
    // This assumption includes the assumption that there are no
    // other instances of StandardHandlerRegistry.
    // TODO make StandardHandlerRegistry a singleton.
    //
    GPTLinitialize();
}



StandardHandlerRegistry::~StandardHandlerRegistry( void )
{
    // We have been using GPTL for data collection, and 
    // since we assume that we are the only GPTL user in the process
    // (see the comment in the ctor), we can gracefully clean up GPTL.
    GPTLfinalize();
}


std::shared_ptr<ITimer>
StandardHandlerRegistry::getTimer(std::string name)
{
    // TODO is there any need for us to retain access to this Timer?
    // TODO do we need to check whether client has already created
    // an object with this name and return that timer?
    return std::make_shared<GPTLTimer>( name );
}


std::shared_ptr<IEventCounter> 
StandardHandlerRegistry::getEventCounter(std::string name)
{
    // TODO is there need for us to retain access to this object?
    // TODO do we need to check whether client has already created
    // an object with this name and return that object?
    return std::make_shared<EventCounter>( name );
}


std::shared_ptr<IHardwareCounter> 
StandardHandlerRegistry::getHardwareCounter( std::string name,
                                            std::vector<HardwareQuantities> quantities)
{
#if READY
#else
    std::cerr << "StandardHandlerRegistry::getHardwareCounter NIY" << std::endl;
#endif // READY
    return std::make_shared<GPTLHardwareCounter>;
}


void
StandardHandlerRegistry::dump(std::ostream& os) const
{
#if READY
#else
    os << "StandardHandlerRegistry::dump NIY" << std::endl;
#endif // READY
}

#endif // RREADY


}//end namespace xolotlPerf

