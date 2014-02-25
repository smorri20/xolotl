#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <GPTLHardwareCounter.h>
#include <vector>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the HardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (HardwareCounter_testSuite)

const std::vector<HardwareQuantities> test_hquan = {L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

BOOST_AUTO_TEST_CASE(checkName) {

	GPTLHardwareCounter tester("test",test_hquan);

	std::cout << "\n" << "GPTLHardwareCounter Message: \n" << "tester.getName() = " << tester.getName() << "\n"
						  << std::endl;

	BOOST_REQUIRE_EQUAL("test", tester.getName());
}

BOOST_AUTO_TEST_CASE(checkInitialValue) {

	GPTLHardwareCounter tester("test",test_hquan);
	std::vector<int> initialValues(test_hquan.size(), 0);

	std::cout << "\n" << "GPTLHardwareCounter Message: \n" << "test_hquan.size() = " << test_hquan.size();

	std::cout << "\n" << "initialValues = ";
	for (std::vector<int>::iterator it = initialValues.begin(); it != initialValues.end(); ++it)
	    std::cout << ' ' << *it;
	std::cout << '\n';

	std::cout << "tester.getValues() = ";
	for (std::vector<int>::iterator it = tester.getValues().begin(); it != tester.getValues().end(); ++it)
	    std::cout << ' ' << *it;
	std::cout << '\n';

	for (unsigned i = 0; i < tester.getValues().size(); i++)
	    BOOST_REQUIRE_EQUAL(initialValues.at(i), tester.getValues().at(i));

}

BOOST_AUTO_TEST_CASE(checkCounting) {

	GPTLHardwareCounter tester("test",test_hquan);

	std::cout << "\n" << "GPTLHardwareCounter Message: \n"
			<< "tester.getValues().size() = " << tester.getValues().size() << std::endl;

//	std::cout << "Values before increment " << " tester.getValues() = ";
//	for (std::vector<int>::iterator it = tester.getValues().begin(); it != tester.getValues().end(); ++it)
//	    std::cout << ' ' << *it;

	tester.increment();

//	std::cout << "\n" << "Values after increment " << " tester.getValues() = ";
//	for (std::vector<int>::iterator it = tester.getValues().begin(); it != tester.getValues().end(); ++it)
//	    std::cout << ' ' << *it;

	std::cout << "\n" << "tester.getValues() = ";
	for (unsigned i = 0; i < tester.getValues().size(); i++)
	{
		//tester.increment();
		std::cout << tester.getValues().at(i) << " ";
	    BOOST_REQUIRE_EQUAL(1, tester.getValues().at(i));
	}

}


BOOST_AUTO_TEST_SUITE_END()
