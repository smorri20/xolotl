#ifndef PLOT_H
#define PLOT_H

// Includes
#include "IPlot.h"
#include "PlottingStyle.h"
#include "DataProvider.h"
#include "LabelProvider.h"
#include <string>

namespace xolotlViz {

/**
 * Plot is the class that realizes the interface IPlot.
 * It is a general class that provides general methods, but to actual plot anything,
 * the user needs to use one of its subclasses.
 */
class Plot: public IPlot {

private:

	/**
	 * Unit of the data that is being plotted.
	 */
	std::string unit;

	/**
	 * Choice of PlottingStyle.
	 */
	PlottingStyle style;

	/**
	 * Label that will appear of the X axis of the plot.
	 */
	std::string axis1Label;

	/**
	 * Label that will appear of the Y axis of the plot.
	 */
	std::string axis2Label;

	/**
	 * Label that will appear of the Z axis of the plot.
	 */
	std::string axis3Label;

	/**
	 * Title of the plot. Can be set by the user by using setTitle().
	 */
	std::string title;

	/**
	 * If it is equal to True, the legend will be displayed.
	 */
	bool enableLegend;

	/**
	 * Data provider used for the plot.
	 */
	DataProvider dataProvider;

	/**
	 * LabelProvider used for the Plot.
	 */
	LabelProvider labelProvider;

public:

	/**
	 * The default constructor
	 */
	Plot();

	/**
	 * The destructor.
	 */
	~Plot();

	/**
	 * Method getting the unit.
	 */
	std::string getUnit();

	/**
	 * Method getting the title.
	 */
	std::string getTitle();

	/**
	 * Method that enables the rendering of the legend.
	 */
	void showLegend(bool legendShow = false);

	/**
	 * Method getting the legend.
	 */
	std::string getLegend();

	/**
	 * Method that will save the plotted plot in a file.
	 */
	void write(std::string fileName);

	/**
	 * Method that allows the user to set his own title.
	 */
	void setTitle(std::string title);

	/**
	 * Method allowing the user to set the PlottingStyle.
	 */
	void setPlottingStyle(PlottingStyle style);

	/**
	 * Sets the unit of the data.
	 */
	void setUnit(std::string unit);

	/**
	 * Sets the data provider used for the plots.
	 */
	void setDataProvider(DataProvider dataProvider);

	/**
	 * Gets the data provider used.
	 */
	DataProvider getDataProvider();

	/**
	 * Sets the label provider used for the plots.
	 */
	void setLabelProvider(LabelProvider labelProvider);

	/**
	 * Gets the label provider used.
	 */
	LabelProvider getLabelProvider();

};

//end class Plot

} /* namespace xolotlViz */

#endif
