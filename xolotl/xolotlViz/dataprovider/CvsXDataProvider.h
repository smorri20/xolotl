#ifndef CVSXDATAPROVIDER_H
#define CVSXDATAPROVIDER_H

// Includes
#include "DataProvider.h"
#include <vector>

namespace xolotlViz {

/**
 * Subclass of DataProvider that will provide the methods to give the value (concentration here)
 * and X data to a ScatterPlot.
 */
class CvsXDataProvider: public DataProvider {

public:

	/**
	 * The default constructor
	 */
	CvsXDataProvider();

	/**
	 * The destructor
	 */
	~CvsXDataProvider();

	/**
	 * Method returning a vector containing the 'x' field of the collection of Point of the DataProvider.
	 */
	virtual std::vector<double> getAxis1Vector() const;

	/**
	 * Method returning a vector containing the 'value' field of the collection of Point of the DataProvider.
	 */
	virtual std::vector<double> getAxis2Vector() const;

};

//end class CvsXDataProvider

} /* namespace xolotlViz */

#endif
