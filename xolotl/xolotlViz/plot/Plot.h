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
	std::string plotUnit;

	/**
	 * Choice of PlottingStyle.
	 */
	PlottingStyle plotStyle;

	/**
	 * Title of the plot. Can be set by the user by using setTitle().
	 */
	std::string plotTitle;

	/**
	 * If it is equal to True, the legend will be displayed.
	 */
	bool enableLegend;

	/**
	 * Data provider used for the plot.
	 */
	DataProvider plotDataProvider;

	/**
	 * LabelProvider used for the Plot.
	 */
	LabelProvider plotLabelProvider;

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
	std::string getUnit() const ;

	/**
	 * Method getting the title.
	 */
	std::string getTitle() const ;

	/**
	 * Method that enables the rendering of the legend.
	 */
	void showLegend(bool legendShow = false);

	/**
	 * Method getting the legend.
	 */
	std::string getLegend() const ;

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
	DataProvider getDataProvider() const ;

	/**
	 * Sets the label provider used for the plots.
	 */
	void setLabelProvider(LabelProvider labelProvider);

	/**
	 * Gets the label provider used.
	 */
	LabelProvider getLabelProvider() const ;

};

//end class Plot

} /* namespace xolotlViz */

#endif
