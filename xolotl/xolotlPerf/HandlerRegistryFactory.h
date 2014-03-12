#ifndef HANDLERREGISTRYFACTORY_H
#define HANDLERREGISTRYFACTORY_H

#include <memory>
#include "IHandlerRegistry.h"

namespace xolotlPerf
{

// Build the desired type of handler registry.
// 'useStdRegistry' indicates whether to build a registry
//      that returns the standard performance objects or the dummy (stub) 
//      objects.
// 
// TODO determine if we need to take an enum instead of a bool,
// if we need to support more than these two types of registries.
//
// Returns true if the handler registry was created successfully.
bool initialize( bool useStdRegistry );

// Access the handler registry.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry( void );

}; // end namespace xolotlPerf

#endif // HANDLERREGISTRYFACTORY_H
