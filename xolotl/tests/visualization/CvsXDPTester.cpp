#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <CvsXDataProvider.h>

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the CvsXDataProvider.
 */
BOOST_AUTO_TEST_SUITE(CvsXDataProvider_testSuite)

/**
 * Method checking the ability to use a name.
 */
BOOST_AUTO_TEST_CASE(checkName) {

	// Create myCvsXDataProvider
	auto myCvsXDataProvider = make_shared<CvsXDataProvider>("myCvsXDataProvider");

	BOOST_REQUIRE_EQUAL("myCvsXDataProvider", myCvsXDataProvider->getName());
}

/**
 * Method checking you can get the axis vectors.
 */
BOOST_AUTO_TEST_CASE(checkGetVector) {

	// Create myCvsXDataProvider
	auto myCvsXDataProvider = make_shared<CvsXDataProvider>("myCvsXDataProvider");

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
	myCvsXDataProvider->setPoints(myPoints);

	// Get the vectors back
	auto axis1Vector = myCvsXDataProvider->getAxis1Vector();
	auto axis2Vector = myCvsXDataProvider->getAxis2Vector();

	// First check the size of the vectors
	BOOST_REQUIRE_EQUAL(axis1Vector.size(), myPoints->size());
	BOOST_REQUIRE_EQUAL(axis2Vector.size(), myPoints->size());

	// Loop on all the points in myPoints
	for (int i = 0; i < myPoints->size(); i++) {
		BOOST_REQUIRE_EQUAL(axis1Vector.at(i), myPoints->at(i).x);
		BOOST_REQUIRE_EQUAL(axis2Vector.at(i), myPoints->at(i).value);
	}
}

BOOST_AUTO_TEST_SUITE_END()
