// Includes
#include "Point.h"

using namespace xolotlViz;

Point::Point() {
}

Point::Point(double v, double tData, double xData, double yData, double zData) {
	value = v;
	t = tData;
	x = xData;
	y = yData;
	z = zData;
}

Point::~Point() {
}
