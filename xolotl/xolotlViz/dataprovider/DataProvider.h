#ifndef DATAPROVIDER_H
#define DATAPROVIDER_H

// Includes
#include "IDataProvider.h"
#include <string>
#include <vector>

namespace xolotlViz {

/**
 *  Realization of the IDataProvider interface. This is a general class with general methods,
 *  to actually get data from the data provider, one needs to use the subclasses.
 */
class DataProvider: public IDataProvider {

protected:

	/**
	 * Collection of data points.
	 */
	std::vector<Point>  data;

public:

	/**
	 * Name of the data contained in the data attribute.
	 */
	std::string dataName;

	/**
	 * Unit of the data contained in the data attribute.
	 */
	std::string dataUnit;

	/**
	 * The default constructor.
	 */
	DataProvider();

	/**
	 * The destructor.
	 */
	~DataProvider();

	/**
	 * Returns a collection of the data points.
	 */
	std::vector<Point> getDataPoints() const;

	/**
	 * Method filling the data collection.
	 */
	void setPoints(std::vector<Point> points);

	/**
	 * Returns the value of the mean of all the data points.
	 */
	double getDataMean() const;

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of quantity that will be plotted on the X axis.
	 * Quantity being x, y, z, t, or value.
	 */
	virtual std::vector<double> getAxis1Vector() const {
		return std::vector<double>();
	}

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of quantity that will be plotted on the Y axis.
	 * Quantity being x, y, z, t, or value.
	 */
	virtual std::vector<double> getAxis2Vector() const {
		return std::vector<double>();
	}

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of quantity that will be plotted on the Z axis.
	 * Quantity being x, y, z, t, or value.
	 */
	virtual std::vector<double> getAxis3Vector() const {
		return std::vector<double>();
	}

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of time steps that will be used for the VideoPlots.
	 */
	virtual std::vector<double> getAxis4Vector() const {
		return std::vector<double>();
	}

};

//end class DataProvider

} /* namespace xolotlViz */
#endif
