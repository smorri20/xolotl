#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "Plot.h"
#include "ScatterPlot.h"
#include "PlottingStyle.h"
#include "DataProvider.h"
#include "LabelProvider.h"

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the Plot class.
 */
BOOST_AUTO_TEST_SUITE(Plot_testSuite)

/**
 * Method checking the ability to choose a PlottingStyle.
 */
BOOST_AUTO_TEST_CASE(checkPlottingStyle) {

	// Create myScatterPlot
	shared_ptr<ScatterPlot> myScatterPlot(
			new ScatterPlot());

	PlottingStyle thePlottingStyle = LINE;

	// Set the PlottingStyle of myScatterPlot
	myScatterPlot->setPlottingStyle(thePlottingStyle);

	// Check it is the right one
	BOOST_REQUIRE_EQUAL(myScatterPlot->getPlottingStyle(), thePlottingStyle);
}

/**
 * Method checking everything related to the data provider.
 */
BOOST_AUTO_TEST_CASE(checkDataProvider) {

	// Create myScatterPlot
	shared_ptr<ScatterPlot> myScatterPlot(
			new ScatterPlot());

	// Create myDataProvider
	shared_ptr<DataProvider> myDataProvider(
			new DataProvider());

	// Create a Point vector
	shared_ptr< vector<xolotlViz::Point> > myPoints(
			new (vector<xolotlViz::Point>));

	// And fill it with some Point
	Point aPoint;
	aPoint.value = 3.; aPoint.t = 1.; aPoint.x = 2.;
	myPoints->push_back(aPoint);
	aPoint.value = 2.; aPoint.t = 3.; aPoint.x = 2.;
	myPoints->push_back(aPoint);
	aPoint.value = 5.; aPoint.t = 6.; aPoint.x = -2.;
	myPoints->push_back(aPoint);
	aPoint.value = -8.; aPoint.t = 8.; aPoint.x = 5.;
	myPoints->push_back(aPoint);
	aPoint.value = -7.; aPoint.t = 7.; aPoint.x = 7.;
	myPoints->push_back(aPoint);

	// Set these points in the myDataProvider
	myDataProvider->setPoints(myPoints);

	// Set myDataProvider in myScatterPlot
	myScatterPlot->setDataProvider(myDataProvider);

	// Methods are missing
    BOOST_FAIL("Methods are missing for checkDataProvider");
}

/**
 * Method checking everything related to the legend.
 */
BOOST_AUTO_TEST_CASE(checkLegend) {
    BOOST_FAIL("checkLegend not implement yet");
}

/**
 * Method checking the writing of the file.
 */
BOOST_AUTO_TEST_CASE(checkWrite) {
    BOOST_FAIL("checkWrite not implement yet");
}

BOOST_AUTO_TEST_SUITE_END()
