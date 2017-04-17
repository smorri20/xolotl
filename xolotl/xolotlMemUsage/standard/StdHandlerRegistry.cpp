#include "mpi.h"
#include <sstream>
#include <iostream>
#include <cassert>
#include <tuple>
#include <cmath>
#include <cstring>
#include "xolotlMemUsage/standard/StdHandlerRegistry.h"
#include "xolotlMemUsage/standard/MemUsageSampler.h"

namespace xolotlMemUsage {

template<typename T, typename V>
void StdHandlerRegistry::CollectAllObjectNames(int myRank,
		const std::map<std::string, std::shared_ptr<T> >& myObjs,
		std::map<std::string, MemUsageObjStatistics<V> >& stats) const {

	// Collect my own object's names.
	std::vector<std::string> myNames;
	CollectMyObjectNames(myObjs, myNames);

    // Figure out how much data we will contribute.
    unsigned int nBytes = 0;
    for(auto currName : myNames) {
        nBytes += (currName.length() + 1);  // 1 extra byte per string for a delimiter
    }

    // Marshal my object's names.
    // We need to be able to unmarshal these also.
    // std::copy into an ostringstream with an ostream_iterator 
    // would be useful, but I don't see how to embed a NUL to 
    // separate the strings.  We could use some other character
    // as a string delimiter, but since we don't restrict what
    // a user can include in their object names, we can't really
    // be sure they didn't use *whatever* character we choose as
    // a delimiter.
	char* myNamesBuf = new char[nBytes];
	char* pName = myNamesBuf;
	for (auto nameIter = myNames.begin(); nameIter != myNames.end();
			++nameIter) {

        size_t currLen = nameIter->length();
		strcpy(pName, nameIter->c_str());
        pName[currLen] = '\0';  // ensure it is NUL-terminated.
		pName += (currLen + 1); // advance to location of next name.
	}
    assert(pName == myNamesBuf + nBytes);

	// Let root know how much space it needs to collect all object names
	unsigned int totalNumBytes = 0;
	MPI_Reduce(&nBytes, &totalNumBytes, 1, MPI_UNSIGNED, MPI_SUM, 0,
			MPI_COMM_WORLD);

	// Provide all names to root.
	// First, provide the amount of data from each process.
	int cwSize;
	MPI_Comm_size(MPI_COMM_WORLD, &cwSize);
	char* allNames = (myRank == 0) ? new char[totalNumBytes] : NULL;
	int* allNameCounts = (myRank == 0) ? new int[cwSize] : NULL;
	int* allNameDispls = (myRank == 0) ? new int[cwSize] : NULL;

	MPI_Gather(&nBytes, 1, MPI_INT, allNameCounts, 1, MPI_INT, 0,
			MPI_COMM_WORLD);

	// Next, root computes the displacements for data from each process.
	if (myRank == 0) {
		allNameDispls[0] = 0;
		for (int i = 1; i < cwSize; ++i) {
			allNameDispls[i] = allNameDispls[i - 1] + allNameCounts[i - 1];
		}
	}

	// Finally, gather all names to the root process.
	MPI_Gatherv(myNamesBuf, nBytes, MPI_CHAR, allNames, allNameCounts,
			allNameDispls, MPI_CHAR, 0, MPI_COMM_WORLD);

	if (myRank == 0) {
		// Process the gathered names to determine the
		// set of all known object names.
		pName = allNames;
		while (pName < (allNames + totalNumBytes)) {
			auto iter = stats.find(pName);
			if (iter == stats.end()) {
				// This is an object  name we have not seen before.
				// Add it to the statistics map.
				stats.insert(
						std::pair<std::string, MemUsageObjStatistics<V> >(pName,
								MemUsageObjStatistics<V>()));
			}

			// Advance to next object name
			pName += (strlen(pName) + 1);
		}
		assert(pName == allNames + totalNumBytes);
	}

	// clean up
	delete[] myNamesBuf;
	delete[] allNames;
	delete[] allNameCounts;
	delete[] allNameDispls;
}

template<typename T>
void StdHandlerRegistry::CollectMyObjectNames(
		const std::map<std::string, std::shared_ptr<T> >& myObjs,
		std::vector<std::string>& objNames) const {
	for (auto oiter = myObjs.begin(); oiter != myObjs.end(); ++oiter) {
		objNames.push_back(oiter->first);
	}
}


template<typename T, typename V>
std::shared_ptr<V> StdHandlerRegistry::GetObjValue(
		const std::map<std::string, std::shared_ptr<T> >& myObjs,
		const std::string& objName) const {

	auto currObjIter = myObjs.find(objName);
	bool found = currObjIter != myObjs.end();

    return found ? std::dynamic_pointer_cast<V>(currObjIter->second->getValue())
                : std::make_shared<V>();
}

template<typename T, typename V>
void StdHandlerRegistry::AggregateStatistics(int myRank,
		const std::map<std::string, std::shared_ptr<T> >& myObjs,
		std::map<std::string, MemUsageObjStatistics<V> >& stats) const {

	// Determine the set of object names known across all processes.
	// Since some processes may define an object that others don't, we
	// have to form the union across all processes.
	// Unfortunately, because the strings are of different lengths,
	// we have a more difficult marshal/unmarshal problem than we'd like.
	CollectAllObjectNames<T, V>(myRank, myObjs, stats);

	// Let all processes know how many statistics we will be collecting.
	int nObjs;
	if (myRank == 0) {
		nObjs = stats.size();
	}
	MPI_Bcast(&nObjs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	assert(nObjs >= 0);

	// Collect and compute statistics for each object.
	auto tsiter = stats.begin();
	for (int idx = 0; idx < nObjs; ++idx) {
		// broadcast the current object's name
		int nameLen = (myRank == 0) ? tsiter->first.length() : -1;
		MPI_Bcast(&nameLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// we can safely cast away const on the tsiter data string because
		// the only process that accesses that string is rank 0,
		// and it only reads the data.
		char* objName =
				(myRank == 0) ?
						const_cast<char*>(tsiter->first.c_str()) :
						new char[nameLen + 1];
		MPI_Bcast(objName, nameLen + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

		// do we know about the current object?
		bool knowObject;
		V myVal;
		std::tie<bool, V>(knowObject, myVal) = GetObjValue<T, V>(myObjs,
				objName);

		// collect count of processes knowing the current object
		unsigned int* pcount =
				(myRank == 0) ? &(tsiter->second.processCount) : NULL;
		int knowObjVal = knowObject ? 1 : 0;
		MPI_Reduce(&knowObjVal, pcount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		// collect min value of current object
		V* pMinVal = (myRank == 0) ? &(tsiter->second.min) : NULL;
		V reduceVal = knowObject ? myVal : T::MaxValue;
		MPI_Reduce(&reduceVal, pMinVal, 1, T::MPIValType, MPI_MIN, 0,
				MPI_COMM_WORLD);

		// collect max value of current object
		V* pMaxVal = (myRank == 0) ? &(tsiter->second.max) : NULL;
		reduceVal = knowObject ? myVal : T::MinValue;
		MPI_Reduce(&reduceVal, pMaxVal, 1, T::MPIValType, MPI_MAX, 0,
				MPI_COMM_WORLD);

		// collect sum of current object's values (for computing avg and stdev)
		double valSum;
		// use the same myVal as for max: actual value if known, 0 otherwise
		double myValAsDouble = (double) reduceVal;
		MPI_Reduce(&myValAsDouble, &valSum, 1, MPI_DOUBLE, MPI_SUM, 0,
				MPI_COMM_WORLD);
		if (myRank == 0) {
			tsiter->second.average = valSum / tsiter->second.processCount;
		}

		// collect sum of squares of current object's values (for stdev)
		double valSquaredSum;
		double myValSquared = myValAsDouble * myValAsDouble;
		MPI_Reduce(&myValSquared, &valSquaredSum, 1, MPI_DOUBLE, MPI_SUM, 0,
				MPI_COMM_WORLD);
		if (myRank == 0) {
			tsiter->second.stdev =
					sqrt(
							(valSquaredSum / tsiter->second.processCount)
									- (tsiter->second.average
											* tsiter->second.average));
		}

		// clean up
		if (myRank != 0) {
			delete[] objName;
		}

		// advance to next object
		if (myRank == 0) {
			++tsiter;
		}
	}
}


template<>
void StdHandlerRegistry::AggregateStatistics<IMemUsageSampler, MemUsageStats>(int myRank,
		const std::map<std::string, std::shared_ptr<IMemUsageSampler> >& myObjs,
		std::map<std::string, MemUsageObjStatistics<MemUsageStats> >& globalStats) const {

	// Determine the set of object names known across all processes.
	// Since some processes may define an object that others don't, we
	// have to form the union across all processes.
	// Unfortunately, because the strings are of different lengths,
	// we have a more difficult marshal/unmarshal problem than we'd like.
	CollectAllObjectNames<IMemUsageSampler, MemUsageStats>(myRank, myObjs, globalStats);

	// Let all processes know how many statistics we will be collecting.
	int nObjs;
	if (myRank == 0) {
		nObjs = globalStats.size();
	}
	MPI_Bcast(&nObjs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	assert(nObjs >= 0);

	// Collect and compute statistics for each object.
	auto tsiter = globalStats.begin();
	for (int idx = 0; idx < nObjs; ++idx) {

		// broadcast the current object's name
		int nameLen = (myRank == 0) ? tsiter->first.length() : -1;
		MPI_Bcast(&nameLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// we can safely cast away const on the tsiter data string because
		// the only process that accesses that string is rank 0,
		// and it only reads the data.
		char* objName =
				(myRank == 0) ?
						const_cast<char*>(tsiter->first.c_str()) :
						new char[nameLen + 1];
		MPI_Bcast(objName, nameLen + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

		// do we know about the current object?
        std::shared_ptr<MemUsageStats> localValue = GetObjValue<IMemUsageSampler, MemUsageStats>(myObjs, objName);
		bool localIsValid = (bool)localValue;

		// collect count of processes knowing about the current object
        uint32_t currProcessCount = 0;
		int knowObjVal = localIsValid ? 1 : 0;
		MPI_Reduce(&knowObjVal, &currProcessCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(myRank == 0) {
            tsiter->second.processCount = currProcessCount;
        }

        // Aggregate values across all ranks.
        tsiter->second.stats.Aggregate(localValue,
                                        myRank,
                                        localIsValid,
                                        currProcessCount);

		// clean up
		if (myRank != 0) {
			delete[] objName;
		}

		// advance to next object
		if (myRank == 0) {
			++tsiter;
		}
	}
}

std::shared_ptr<IHandlerRegistry::GlobalData>
StdHandlerRegistry::collectData(void) const {

    auto ret = std::make_shared<GlobalMemUsageStats>();

	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	// Aggregate data from all processes.
    // Memory usage.
    AggregateStatistics<IMemUsageSampler, MemUsageStats>(myRank,
        allSamplers, ret->memStats);

    return ret;
}

void StdHandlerRegistry::reportData(std::ostream& os,
                        std::shared_ptr<IHandlerRegistry::GlobalData> data) const {

    auto myData = std::dynamic_pointer_cast<GlobalMemUsageStats>(data);

    os << "\nMemoryUsage:\n";
    for (auto currMemUsageStats : myData->memStats) {
        os << "name: " << currMemUsageStats.first << '\n';
        currMemUsageStats.second.outputTo(os);
        os << '\n';
    }
}

} // namespace xolotlMemUsage

