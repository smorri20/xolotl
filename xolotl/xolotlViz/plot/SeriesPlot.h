#ifndef SERIESPLOT_H
#define SERIESPLOT_H

// Includes
#include "Plot.h"
#include <vector>

namespace xolotlViz {

/**
 * Plot the different data values as a function of one dimension. Available PlottingStyle are POINTS or LINE.
 * It can be associated to CvsXDataProvider.
 */
class SeriesPlot: public Plot {

private:

	/**
	 * Container of data providers used for the plot.
	 */
	std::shared_ptr< std::vector< std::shared_ptr<DataProvider> > > plotDataProviders;

public:

	/**
	 * The default constructor
	 */
	SeriesPlot();

	/**
	 * The destructor
	 */
	~SeriesPlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	void render();

	/**
	 * Method adding one data provider to the vector plotDataProviders
	 */
	void addDataProvider(std::shared_ptr<DataProvider> dataProvider);

	/**
	 * Method getting the i-th data provider
	 */
	std::shared_ptr<DataProvider> getDataProvider(int i);

	/**
	 * Method getting the total number of data providers
	 */
	int getDataProviderNumber();

	/**
	 * Method getting the maximum value taken by the data
	 */
	double getMaxValue();

};

//end class SeriesPlot

} /* namespace xolotlViz */

#endif
