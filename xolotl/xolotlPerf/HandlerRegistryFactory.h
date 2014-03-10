#ifndef PERFFACTORY_H
#define PERFFACTORY_H

#include <memory>
#include "IHandlerRegistry.h"

namespace xolotlPerf
{

// Based on command line arguments, build the desired type of handler registry.
void initialize( int argc, char* argv[] );

// Access the handler registry.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry( void );

}; // end namespace xolotlPerf

#endif // PERFFACTORY_H
