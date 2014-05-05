#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyPlot.h>
#include <PlottingStyle.h>
#include <DataProvider.h>

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the DummyPlot class.
 */
BOOST_AUTO_TEST_SUITE(DummyPlot_testSuite)

/**
 * Method checking the non-ability to use a name.
 */
BOOST_AUTO_TEST_CASE(checkName) {

	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>("myDummyPlot");

	BOOST_REQUIRE_EQUAL("unused", myDummyPlot->getName());
}

/**
 * Method checking the non-ability to choose a PlottingStyle.
 */
BOOST_AUTO_TEST_CASE(checkPlottingStyle) {

	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>("myDummyPlot");

	PlottingStyle thePlottingStyle = LINE;

	// Set the PlottingStyle of myDummyPlot
	myDummyPlot->setPlottingStyle(thePlottingStyle);

	// Check it is the right one
	BOOST_REQUIRE_EQUAL(myDummyPlot->getPlottingStyle(), PlottingStyle());
}

/**
 * Method checking everything related to the data provider.
 */
BOOST_AUTO_TEST_CASE(checkDataProvider) {

	// Create myDummyPlot
	auto myDummyPlot = make_shared<DummyPlot>("myDummyPlot");

	// Create myDataProvider
	auto myDataProvider = make_shared<DataProvider>("myDataProvider");

	// Create a Point vector
	auto myPoints = make_shared< vector <Point> >();

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

	// Set myDataProvider in myDummyPlot
	myDummyPlot->setDataProvider(myDataProvider);

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
