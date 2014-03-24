#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DataProvider.h>

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the DataProvider.
 */
BOOST_AUTO_TEST_SUITE(DataProvider_testSuite)

/**
 * Method checking you can set and get the unit of the data.
 */
BOOST_AUTO_TEST_CASE(checkUnit) {
	// Create myDataProvider
	shared_ptr<DataProvider> myDataProvider(
			new DataProvider());

	string theUnit = "mol per cubic meter";

	myDataProvider->dataUnit = theUnit;

	BOOST_REQUIRE_EQUAL(myDataProvider->dataUnit, theUnit);
}

/**
 * Method checking you can set and get the name of the data.
 */
BOOST_AUTO_TEST_CASE(checkName) {
	// Create myDataProvider
	shared_ptr<DataProvider> myDataProvider(
			new DataProvider());

	string theName = "Concentration";

	myDataProvider->dataName = theName;

	BOOST_REQUIRE_EQUAL(myDataProvider->dataName, theName);
}

/**
 * Method checking you can add points to the data, get the data, and getDataMean().
 */
BOOST_AUTO_TEST_CASE(checkData) {
	// Create myDataProvider
	shared_ptr<DataProvider> myDataProvider(
			new DataProvider());

	// Create a Point vector
	vector<xolotlViz::Point> myPoints;

	myPoints.push_back(Point(3., 1., 2.));
	myPoints.push_back(Point(2., 3., 2., 4.));
	myPoints.push_back(Point(5., 6., -2., 7., -1.));
	myPoints.push_back(Point(-8., 8., 5.));
	myPoints.push_back(Point(7., 7., 7.));

	// Set these points in the myDataProvider
	myDataProvider->setPoints(myPoints);

	auto dataPoints = myDataProvider->getDataPoints();

	// First check the size of the vector
	BOOST_REQUIRE_EQUAL(dataPoints.size(), myPoints.size());

	// Loop on all the points in dataPoints
	for (int i = 0; i < dataPoints.size(); i++) {
		BOOST_REQUIRE_EQUAL(dataPoints.at(i).value, myPoints.at(i).value);
		BOOST_REQUIRE_EQUAL(dataPoints.at(i).t, myPoints.at(i).t);
		BOOST_REQUIRE_EQUAL(dataPoints.at(i).x, myPoints.at(i).x);
		BOOST_REQUIRE_EQUAL(dataPoints.at(i).y, myPoints.at(i).y);
		BOOST_REQUIRE_EQUAL(dataPoints.at(i).z, myPoints.at(i).z);
	}

	// Get the mean value of the data
	auto mean = myDataProvider->getDataMean();

	BOOST_REQUIRE_EQUAL(mean, 1.8);
}

BOOST_AUTO_TEST_SUITE_END()
