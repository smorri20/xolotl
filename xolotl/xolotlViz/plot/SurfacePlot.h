#ifndef SURFACEPLOT_H
#define SURFACEPLOT_H

// Includes
#include "Plot.h"
#include "Point.h"
#include <vector>

namespace xolotlViz {

/**
 * Plot the data value as a function of two different dimensions.
 * The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE.
 * It can be associated to QvsXYDataProvider, QvsXZDataProvider, or QvsYZDataProvider.
 */
class SurfacePlot: public Plot {

private:

	/**
	 * Vector containing three fields: one for the value and two for the directions as a function
	 * of which the value is plotted.
	 */
	Vector data;

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
	void render(bool isCumulative = false);

	/**
	 * Method returning the data points that are stored in the data vector.
	 */
	Point getPoints();

	/**
	 * Method getting the X axis label with the help of the label provider.
	 */
	std::string getAxis1Label();

	/**
	 * Method getting the Y axis label with the help of the label provider.
	 */
	std::string getAxis2Label();

	/**
	 * Method getting the Z axis label with the help of the label provider.
	 */
	std::string getAxis3Label();

};

//end class SurfacePlot

} /* namespace xolotlViz */

#endif
