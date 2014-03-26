#ifndef SCATTERPLOT_H
#define SCATTERPLOT_H

// Includes
#include "Plot.h"
#include <vector>

namespace xolotlViz {

/**
 * Plot the data value as a function of one dimension. Available PlottingStyle are POINTS or LINE.
 * It can be associated to QvsXDataProvider, QvsYDataProvider, QvsZDataProvider, or QvsTimeDataProvider.
 */
class ScatterPlot: public Plot {

public:

	/**
	 * The default constructor
	 */
	ScatterPlot();

	/**
	 * The destructor
	 */
	~ScatterPlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	void render();

};

//end class ScatterPlot

} /* namespace xolotlViz */

#endif
