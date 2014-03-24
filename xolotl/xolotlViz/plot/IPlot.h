#ifndef IPLOT_H
#define IPLOT_H

// Includes
#include <string>

namespace xolotlViz {

/**
 * IPlot describe the structure needed to be able to plot data provided by IDataProvider.
 * The user interacts with it through different method where he/she could set the data to plot,
 * title, legend, plotting style, etc.
 */
class IPlot {

public:

	/**
	 * Method that enables the rendering of the legend.
	 */
	virtual void showLegend(bool legendShow = false) = 0;

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	virtual void render(bool isCumulative = false) = 0;

	/**
	 * Method that will save the plotted plot in a file.
	 */
	virtual void write(std::string fileName) = 0;

	/**
	 * Method that allows the user to set his own title.
	 */
	virtual void setTitle(std::string title) = 0;

	/**
	 * Method allowing the user to set the PlottingStyle.
	 */
	virtual void setPlottingStyle(PlottingStyle style) = 0;

	/**
	 * Method allowing the user to set the unit of the value plotted.
	 */
	virtual void setUnit(std::string unit) = 0;

	/**
	 * Sets the data provider used for the plots.
	 */
	virtual void setDataProvider(DataProvider dataProvider) = 0;

	/**
	 * Gets the data provider used.
	 */
	virtual DataProvider getDataProvider() = 0;

	/**
	 * Sets the label provider used for the plots.
	 */
	virtual void setLabelProvider(LabelProvider labelProvider) = 0;

	/**
	 * Gets the label provider used.
	 */
	virtual LabelProvider getLabelProvider() = 0;

	/**
	 * Method getting the unit from the data provider and setting the attribute unit with this value.
	 */
	virtual std::string getUnit() = 0;

	/**
	 * Method getting the title with the help of the data provider and the label provider and
	 * setting the title attribute.
	 */
	virtual std::string getTitle() = 0;

	/**
	 * Method defining the legend with the help of the data provider and the label provider.
	 */
	virtual std::string getLegend() = 0;

};

//end class IPlot

} /* namespace xolotlViz */

#endif
