#ifndef IVIZHANDLERREGISTRY_H
#define IVIZHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "IPlot.h"
#include "PlotType.h"


namespace xolotlViz {

/**
 * Factory for building plots, dataprovider, and labelprovider.
 */
class IVizHandlerRegistry {

public:

    /// Possible types of visualization handler registries.
    enum RegistryType
    {
        dummy,      //< Use stub classes that do nothing
        std,        //< Use the best available API.
        png,        //< Use raster based rendering and output to raster PNG file
        eps         //< Use vector based rendering and output to vector EPS file
    };

	/**
	 * The destructor
	 */
	virtual ~IVizHandlerRegistry(){}

	/**
	 * This operation returns the IPlot specified by the parameter.
	 */
	virtual std::shared_ptr<IPlot> getPlot(std::string name, PlotType type) = 0;

}; //end class IVizHandlerRegistry

} //end namespace xolotlViz

#endif
