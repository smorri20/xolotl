#ifndef SURFACEPLOT_H
#define SURFACEPLOT_H

// Includes
#include "Plot.h"
#include <vector>

namespace xolotlViz {

/**
 * Plot the data value as a function of two different dimensions.
 * The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE.
 * It can be associated to QvsXYDataProvider, QvsXZDataProvider, or QvsYZDataProvider.
 */
class SurfacePlot: public Plot {

public:

	/**
	 * The default constructor
	 */
	SurfacePlot();

	/**
	 * The destructor
	 */
	~SurfacePlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	void render();

};

//end class SurfacePlot

} /* namespace xolotlViz */

#endif
