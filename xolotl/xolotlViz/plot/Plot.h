#ifndef PLOT_H
#define PLOT_H

// Includes
#include "IPlot.h"
#include "PlottingStyle.h"
#include "../dataprovider/DataProvider.h"
#include "../labelprovider/LabelProvider.h"

namespace xolotlViz {

/**
 * Plot is the class that realizes the interface IPlot.
 * It is a general class that provides general methods, but to actual plot anything,
 * the user needs to use one of its subclasses.
 */
class Plot: public IPlot {

protected:

	/**
	 * Choice of PlottingStyle.
	 */
	PlottingStyle plotStyle;

	/**
	 * If it is equal to True, the legend will be displayed.
	 */
	bool enableLegend = false;

	/**
	 * If it is equal to True, a log scale will be used (for 1D plot for now).
	 */
	bool enableLogScale = false;

	/**
	 * Data provider used for the plot.
	 */
	std::shared_ptr<DataProvider> plotDataProvider;

public:

	/**
	 * LabelProvider used for the Plot.
	 */
	std::shared_ptr<LabelProvider> plotLabelProvider;

	/**
	 * The default constructor
	 */
	Plot();

	/**
	 * The destructor.
	 */
	~Plot();

	/**
	 * Method that will save the plotted plot in a file.
	 * \see IPlot.h
	 */
	void write(std::string fileName);

	/**
	 * Method allowing the user to set the PlottingStyle.
	 * \see IPlot.h
	 */
	void setPlottingStyle(PlottingStyle style);

	/**
	 * Method getting the PlottingStyle.
	 * \see IPlot.h
	 */
	PlottingStyle getPlottingStyle();

	/**
	 * Sets the data provider used for the plots.
	 * \see IPlot.h
	 */
	void setDataProvider(std::shared_ptr<DataProvider> dataProvider);

	/**
	 * Gets the data provider used.
	 * \see IPlot.h
	 */
	std::shared_ptr<DataProvider> getDataProvider() const ;

	/**
	 * Sets the label provider used for the plots.
	 * \see IPlot.h
	 */
	void setLabelProvider(std::shared_ptr<LabelProvider> labelProvider);

	/**
	 * Gets the label provider used.
	 * \see IPlot.h
	 */
	std::shared_ptr<LabelProvider> getLabelProvider() const ;

	/**
	 * Method that enables the rendering of the legend.
	 * \see IPlot.h
	 */
	void showLegend(bool legendShow = true);

	/**
	 * Method getting the legend.
	 * \see IPlot.h
	 */
	std::string getLegend() const ;

	/**
	 * Method that enables the log scale.
	 * \see IPlot.h
	 */
	void setLogScale(bool logScale = true);

};

//end class Plot

} /* namespace xolotlViz */

#endif
