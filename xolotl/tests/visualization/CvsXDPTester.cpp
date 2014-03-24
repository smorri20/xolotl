#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <CvsXDataProvider.h>

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the DataProvider.
 */
BOOST_AUTO_TEST_SUITE(DataProvider_testSuite)

/**
 * Method checking you can get the axis vectors.
 */
BOOST_AUTO_TEST_CASE(checkGetVector) {
	// Create myDataProvider
	shared_ptr<CvsXDataProvider> myCvsXDataProvider(
			new CvsXDataProvider());

	// Create a Point vector
	vector<xolotlViz::Point> myPoints;

	myPoints.push_back(Point(3., 1., 2.));
	myPoints.push_back(Point(2., 3., 2., 4.));
	myPoints.push_back(Point(5., 6., -2., 7., -1.));
	myPoints.push_back(Point(-8., 8., 5.));
	myPoints.push_back(Point(7., 7., 7.));

	// Set these points in the myDataProvider
	myCvsXDataProvider->setPoints(myPoints);

	// Get the vectors back
	auto axis1Vector = myCvsXDataProvider->getAxis1Vector();
	auto axis2Vector = myCvsXDataProvider->getAxis2Vector();

	// First check the size of the vectors
	BOOST_REQUIRE_EQUAL(axis1Vector.size(), myPoints.size());
	BOOST_REQUIRE_EQUAL(axis2Vector.size(), myPoints.size());

	// Loop on all the points in myPoints
	for (int i = 0; i < myPoints.size(); i++) {
		BOOST_REQUIRE_EQUAL(axis1Vector.at(i), myPoints.at(i).x);
		BOOST_REQUIRE_EQUAL(axis2Vector.at(i), myPoints.at(i).value);
	}
}

BOOST_AUTO_TEST_SUITE_END()
