#include <cassert>
#include "StandardHandlerRegistry.h"
#include "gptl.h"


namespace xolotlPerf
{


StandardHandlerRegistry::StandardHandlerRegistry( std::vector<HardwareQuantities> hwq )
{
    // We use GPTL for data collection, so we must make sure 
    // that library has been initialized.
    //
    // Note: We assume that no other part of the code is using GPTL.
    // This assumption includes the assumption that there are no
    // other instances of StandardHandlerRegistry.
    //

    // Indicate to GPTL any hardware counters it should be monitoring.
    for( auto iter = hwq.begin(); iter != hwq.end(); iter++ )
    {
        HardwareQuantities currHwq = *iter;

        // Convert the current HardwareQuantity into PAPI notation,
        // since GPTL doesn't understand our enum values.
        // The value should be found in the map.
        HardwareQuantityInfoMap::const_iterator mapiter = hwqInfoMap.find( currHwq );
        assert( mapiter != hwqInfoMap.end() );
        const HardwareQuantityInfo& currHwqInfo = mapiter->second;

        int gret = GPTLsetoption( currHwqInfo.papiID, 1 );
        if( gret < 0 )
        {
            std::cerr << "Warning: unable to monitor requested hardware counter " 
                << currHwqInfo.name
                << " (PAPI name: " << currHwqInfo.papiName
                << ')'
                << std::endl;
        }
    }

    // Initialize the GPTL library.
    // GPTLsetoption(GPTLverbose, 1);   // useful for debugging
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
    // TODO is there need for us to retain access to this object?
    // TODO do we need to check whether this client has already
    // created an object with this name and return that object?
    return std::make_shared<GPTLHardwareCounter>( name, quantities );
}


void
StandardHandlerRegistry::dump(std::ostream& os) const
{
#if READY
#else
    os << "StandardHandlerRegistry::dump NIY" << std::endl;
#endif // READY
}


}//end namespace xolotlPerf

