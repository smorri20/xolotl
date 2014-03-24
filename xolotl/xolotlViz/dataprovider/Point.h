#ifndef POINT_H
#define POINT_H

// Includes
#include <string>

namespace xolotlViz {

/**
 * Class describing the structure of data points.
 * The attributes are the three spatial dimensions, the time, and the value of the quantity
 * under consideration at this position.
 */
class Point {


public:

	/**
	 * The time step.
	 */
	double t;

	/**
	 * The X position on the grid.
	 */
	double x;

	/**
	 * The Y position on the grid.
	 */
	double y;

	/**
	 * The Z position on the grid.
	 */
	double z;

	/**
	 * Value of the quantity of interest at the time step t and position on the grid x,y,z.
	 */
	double value;

	/**
	 * The default constructor
	 */
	Point();

	/**
	 * Constructor initializing the values of the point.
	 * @param v will fill the value attribute of the Point.
	 * @param tData will fill the t attribute of the Point.
	 * @param xData will fill the x attribute of the Point.
	 * @param yData will fill the y attribute of the Point.
	 * @param zData will fill the z attribute of the Point.
	 */
	Point(double v, double tData, double xData, double yData = 0., double zData = 0.);

	/**
	 * The destructor
	 */
	~Point();

};

//end class Point

} /* namespace xolotlViz */

#endif
