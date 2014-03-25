#ifndef VIDEOPLOT_H
#define VIDEOPLOT_H

// Includes
#include "Plot.h"
#include "Point.h"
#include <vector>

namespace xolotlViz {

/**
 * Plot the data value as a function of two different spatial dimensions for each time step
 * and change the time step with time to have a video-like rendering, each frame being a SurfacePlot.
 * The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE.
 * It can be associated to QvsXYTimeDataProvider, QvsXZTimeDataProvider, or QvsYZTimeDataProvider.
 */
class VideoPlot: public Plot {

private:

	/**
	 * Vector containing four fields: one for the value, two for the spatial directions
	 * as a function of which the value is plotted, and the last one for the time step.
	 */
	Vector data;

	/**
	 * Number of frames shown per second.
	 */
	double frameRate;

public:

	/**
	 * The default constructor
	 */
	VideoPlot();

	/**
	 * The destructor
	 */
	~VideoPlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	void render();

	/**
	 * Sets the frame rate for VideoPlot.
	 */
	void setFrameRate(double frameRate);

	/**
	 * Get the frame rate from VideoPlot.
	 */
	double getFrameRate();

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

	/**
	 * Method getting the time axis label with the help of the label provider.
	 */
	std::string getAxis4Label();

};

//end class VideoPlot

} /* namespace xolotlViz */

#endif
