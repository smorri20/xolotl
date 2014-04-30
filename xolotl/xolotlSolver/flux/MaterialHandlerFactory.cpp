#include "MaterialHandlerFactory.h"
#include "WFitFluxHandler.h"
#include "FeFitFluxHandler.h"
#include <XolotlOptions.h>
#include <iostream>

namespace xolotlSolver
{

static std::shared_ptr<IFluxHandler> theFluxHandler;

// Create the desired type of handler registry.
bool initializeMaterial( bool useWRegistry )
{
    bool ret = true;

    if( useWRegistry )
    {
        // we are to use a tungsten, W, flux handler
        theFluxHandler = std::make_shared<WFitFluxHandler>( );
    }
    else
    {
        // use a iron, Fe, FitFluxHandler for this run
        theFluxHandler = std::make_shared<FeFitFluxHandler>();
    }

    return ret;
}

// Provide access to our handler registry.
std::shared_ptr<IFluxHandler> getMaterialHandler( void )
{
    if( !theFluxHandler )
    {
        // We have not yet been initialized.
        // Issue a warning and use the default tungsten registry.
        std::cerr << "Warning: xolotlSolver material handler requested, but "
        		"library has not been initialized.  Defaulting to using "
        		"tungsten (W) handlers" << std::endl;

        xolotlSolver::initializeMaterial( true );
    }
    return theFluxHandler;
}

} // end namespace xolotlSolver

