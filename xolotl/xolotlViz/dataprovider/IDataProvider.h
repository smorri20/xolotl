#ifndef IDATAPROVIDER_H
#define IDATAPROVIDER_H

// Includes
#include "Point.h"
#include <vector>
#include <string>

namespace xolotlViz {

/**
 * IDataProvider describes the structure needed as a link between Xololt outputs
 * and quantities that the user wants to plot.
 */
class IDataProvider {

public:

	/**
	 * Returns a collection of the data points.
	 */
	virtual std::vector<Point> getDataPoints() const = 0;

	/**
	 * Returns the value of the mean of all the data points.
	 */
	virtual double getDataMean() const = 0;

	/**
	 * Method filling the data collection.
	 */
	virtual void setPoints(std::vector<Point> Points) = 0;

	/**
	 * Method returning the vector of quantity that will be plotted on the X axis.
	 * Quantity being x, y, z, t, or value.
	 */
	virtual std::vector<double> getAxis1Vector() const = 0;

	/**
	 * Method returning the vector of quantity that will be plotted on the Y axis.
	 * Quantity being x, y, z, t, or value.
	 */
	virtual std::vector<double> getAxis2Vector() const = 0;

	/**
	 * Method returning the vector of quantity that will be plotted on the Z axis.
	 * Quantity being x, y, z, t, or value.
	 */
	virtual std::vector<double> getAxis3Vector() const = 0;

	/**
	 * Method returning the vector of time steps that will be used for the VideoPlots.
	 */
	virtual std::vector<double> getAxis4Vector() const = 0;

};

//end class IDataProvider

} /* namespace xolotlViz */
#endif
