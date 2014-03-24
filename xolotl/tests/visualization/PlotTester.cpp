#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "Plot.h"

using namespace std;
using namespace xolotlViz;

/**
 * This suite is responsible for testing the Plot class.
 */
BOOST_AUTO_TEST_SUITE(Plot_testSuite)

/**
 * Method checking getting and setting the unit.
 */
BOOST_AUTO_TEST_CASE(checkUnit) {
}

/**
 * Method checking the ability to choose a PlottingStyle.
 */
BOOST_AUTO_TEST_CASE(checkPlottingStyle) {
}

/**
 * Method checking everything related to the label provider.
 */
BOOST_AUTO_TEST_CASE(checkLabelProvider) {
}

/**
 * Method checking the ability to change the title of the plot.
 */
BOOST_AUTO_TEST_CASE(checkTitle) {
}

/**
 * Method checking everything related to the legend.
 */
BOOST_AUTO_TEST_CASE(checkLegend) {
}

/**
 * Method checking the writing of the file.
 */
BOOST_AUTO_TEST_CASE(checkWrite) {
}

/**
 * Method checking everything related to the data provider.
 */
BOOST_AUTO_TEST_CASE(checkDataProvider) {
}

BOOST_AUTO_TEST_SUITE_END()
