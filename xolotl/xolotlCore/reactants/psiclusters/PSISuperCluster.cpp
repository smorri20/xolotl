// Includes
#include "xolotlCore/reactants/reactantsConfig.h"
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <MathUtils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

#if READY
static std::shared_ptr<xolotlPerf::IEventCounter> resultFrom_addsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> resultFrom_callsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> resultFrom_baddsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> resultFrom_bcallsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> partInProd_addsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> partInProd_callsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> partInProd_baddsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> partInProd_bcallsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> partInDiss_addsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> partInDiss_callsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> partInDiss_baddsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> partInDiss_bcallsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> emitFrom_addsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> emitFrom_callsCounter;

static std::shared_ptr<xolotlPerf::IEventCounter> emitFrom_baddsCounter;
static std::shared_ptr<xolotlPerf::IEventCounter> emitFrom_bcallsCounter;
#endif // READY

/**
 * The helium momentum partials.
 */
std::vector<double> heMomentumPartials;

/**
 * The vacancy momentum partials.
 */
std::vector<double> vMomentumPartials;

PSISuperCluster::PSISuperCluster(double _numHe, double _numV, int _nTot,
		int heWidth, int vWidth,
        IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(_numHe, _numV)),
        numHe(_numHe), numV(_numV), nTot(_nTot),
        heBounds(0, 0),
        vBounds(0, 0),
        l0(0.0), l1He(0.0), l1V(0.0), dispersionHe(
				0.0), dispersionV(0.0), heMomentumFlux(0.0), vMomentumFlux(0.0) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = (int) (numHe + numV);

	// Update the composition map
	composition[toCompIdx(Species::He)] = (int) numHe;
	composition[toCompIdx(Species::V)] = (int) numV;

	// Set the width
	sectionHeWidth = heWidth;
	sectionVWidth = vWidth;

	// Set the formation energy
	formationEnergy = 0.0; // It is set to 0.0 because we do not want the super clusters to undergo dissociation

	// Set the diffusion factor and the migration energy
	migrationEnergy = std::numeric_limits<double>::infinity();
	diffusionFactor = 0.0;

	// Set the typename appropriately
	type = ReactantType::PSISuper;

	return;
}


void PSISuperCluster::resultFrom(ProductionReaction& reaction,
                        int a, int b, int c, int d) {

#if READY
    if (not resultFrom_callsCounter) {
        resultFrom_callsCounter = handlerRegistry->getEventCounter("PSISuper_resultFrom_calls");
    }
    resultFrom_callsCounter->increment();
#endif // READY

	// Check if we already know about the reaction.
    auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
    auto it = effReactingList.find(rkey);
	if (it == effReactingList.end()) {

#if READY
        if (not resultFrom_addsCounter) {
            resultFrom_addsCounter = handlerRegistry->getEventCounter("PSISuper_resultFrom_adds");
        }
        resultFrom_addsCounter->increment();
#endif // READY

		// We did not already know about this reaction.
        // Add info about production to our list.
        auto eret = effReactingList.emplace(std::piecewise_construct,
                            std::forward_as_tuple(rkey),
                            std::forward_as_tuple(
                                reaction,
                                static_cast<PSICluster&>(reaction.first),
                                static_cast<PSICluster&>(reaction.second)));
        // Since we already checked and didn't know about the reaction,
        // we had better have added it with our emplace() call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effReactingList.end());
    auto& prodPair = it->second;    

    // NB: prodPair's reactants are same as reaction.
    // So use prodPair only from here on.
    // TODO any way to enforce this?

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
			secondVDistance = 0.0;
	if (prodPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.first);
		firstHeDistance = super.getHeDistance(c);
		firstVDistance = super.getVDistance(d);
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(prodPair.second);
		secondHeDistance = super.getHeDistance(c);
		secondVDistance = super.getVDistance(d);
	}
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vFactor = (double) (b - numV) / dispersionV;
	// First is A, second is B, in A + B -> this
	prodPair.a[0][0][0] += 1.0;
	prodPair.a[0][0][1] += heFactor;
	prodPair.a[0][0][2] += vFactor;
	prodPair.a[1][0][0] += firstHeDistance;
	prodPair.a[1][0][1] += firstHeDistance * heFactor;
	prodPair.a[1][0][2] += firstHeDistance * vFactor;
	prodPair.a[2][0][0] += firstVDistance;
	prodPair.a[2][0][1] += firstVDistance * heFactor;
	prodPair.a[2][0][2] += firstVDistance * vFactor;
	prodPair.a[0][1][0] += secondHeDistance;
	prodPair.a[0][1][1] += secondHeDistance * heFactor;
	prodPair.a[0][1][2] += secondHeDistance * vFactor;
	prodPair.a[0][2][0] += secondVDistance;
	prodPair.a[0][2][1] += secondVDistance * heFactor;
	prodPair.a[0][2][2] += secondVDistance * vFactor;
	prodPair.a[1][1][0] += firstHeDistance * secondHeDistance;
	prodPair.a[1][1][1] += firstHeDistance * secondHeDistance * heFactor;
	prodPair.a[1][1][2] += firstHeDistance * secondHeDistance * vFactor;
	prodPair.a[1][2][0] += firstHeDistance * secondVDistance;
	prodPair.a[1][2][1] += firstHeDistance * secondVDistance * heFactor;
	prodPair.a[1][2][2] += firstHeDistance * secondVDistance * vFactor;
	prodPair.a[2][1][0] += firstVDistance * secondHeDistance;
	prodPair.a[2][1][1] += firstVDistance * secondHeDistance * heFactor;
	prodPair.a[2][1][2] += firstVDistance * secondHeDistance * vFactor;
	prodPair.a[2][2][0] += firstVDistance * secondVDistance;
	prodPair.a[2][2][1] += firstVDistance * secondVDistance * heFactor;
	prodPair.a[2][2][2] += firstVDistance * secondVDistance * vFactor;

	return;
}

void PSISuperCluster::resultFrom(ProductionReaction& reaction,
        const std::vector<PendingProductionReactionInfo>& prInfos) {

#if READY
    if (not resultFrom_bcallsCounter) {
        resultFrom_bcallsCounter = handlerRegistry->getEventCounter("PSISuper_resultFrom_bcalls");
    }
    resultFrom_bcallsCounter->increment();
#endif // READY

	// Check if we already know about the reaction.
    auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
    auto it = effReactingList.find(rkey);
	if (it == effReactingList.end()) {

#if READY
        if (not resultFrom_baddsCounter) {
            resultFrom_baddsCounter = handlerRegistry->getEventCounter("PSISuper_resultFrom_badds");
        }
        resultFrom_baddsCounter->increment();
#endif // READY

		// We did not already know about this reaction.
        // Add info about production to our list.
        auto eret = effReactingList.emplace(std::piecewise_construct,
                            std::forward_as_tuple(rkey),
                            std::forward_as_tuple(
                                reaction,
                                static_cast<PSICluster&>(reaction.first),
                                static_cast<PSICluster&>(reaction.second)));
        // Since we already checked and didn't know about the reaction,
        // we had better have added it with our emplace() call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effReactingList.end());
    auto& prodPair = it->second;    

    // NB: prodPair's reactants are same as reaction.
    // So use prodPair only from here on.
    // TODO any way to enforce this?

	// Update the coefficients
#if defined(USE_SORTED_COEFF_SUMS)

    // Update coefficients in a data-parallel way.
    // This requires more memory than the original approach since we 
    // store all of the values computed from the prInfos so that we can 
    // sort them before summing them for each (i,j,k) index in prodPair.a.

    // Allocate space for all terms.
    // We do one large flat allocation, then drop "pointers" (iterators)
    // into it at appropriate offsets for each (i,j,k) triple.
    // This is more efficient than doing 27 memory allocations each
    // time we are called, and lets us still use the iterator-based C++
    // standard library functions like std::sort.
    std::vector<double> coeffTermsFlat(prInfos.size()*27);
    Array3D<std::pair<std::vector<double>::iterator, std::vector<double>::iterator>, 3, 3, 3> coeffTerms;
    auto currIter = coeffTermsFlat.begin();
    for(auto i = 0; i < 3; ++i)
    {
        for(auto j = 0; j < 3; ++j)
        {
            for(auto k = 0; k < 3; ++k)
            {
                auto nextIter = currIter + prInfos.size();
                coeffTerms[i][j][k] = std::make_pair(currIter, nextIter);
                currIter = nextIter;
            }
        }
    }

    // Compute the values associated with each PRI.
    auto currIdx = 0;
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&prodPair,&coeffTerms,&currIdx](const PendingProductionReactionInfo& currPRI) {

        // Use names corresponding to those in single version.
        int a = currPRI.numHe;
        int b = currPRI.numV;
        int c = currPRI.i;
        int d = currPRI.j;

        double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
                secondVDistance = 0.0;
        if (prodPair.first.getType() == ReactantType::PSISuper) {
            auto const& super = static_cast<PSICluster const&>(prodPair.first);
            firstHeDistance = super.getHeDistance(c);
            firstVDistance = super.getVDistance(d);
        }
        if (prodPair.second.getType() == ReactantType::PSISuper) {
            auto const& super = static_cast<PSICluster const&>(prodPair.second);
            secondHeDistance = super.getHeDistance(c);
            secondVDistance = super.getVDistance(d);
        }
        double heFactor = (double) (a - numHe) / dispersionHe;
        double vFactor = (double) (b - numV) / dispersionV;
        // First is A, second is B, in A + B -> this
        *(coeffTerms[0][0][0].first + currIdx) = 1.0;
        *(coeffTerms[0][0][1].first + currIdx) = heFactor;
        *(coeffTerms[0][0][2].first + currIdx) = vFactor;
        *(coeffTerms[1][0][0].first + currIdx) = firstHeDistance;
        *(coeffTerms[1][0][1].first + currIdx) = firstHeDistance * heFactor;
        *(coeffTerms[1][0][2].first + currIdx) = firstHeDistance * vFactor;
        *(coeffTerms[2][0][0].first + currIdx) = firstVDistance;
        *(coeffTerms[2][0][1].first + currIdx) = firstVDistance * heFactor;
        *(coeffTerms[2][0][2].first + currIdx) = firstVDistance * vFactor;
        *(coeffTerms[0][1][0].first + currIdx) = secondHeDistance;
        *(coeffTerms[0][1][1].first + currIdx) = secondHeDistance * heFactor;
        *(coeffTerms[0][1][2].first + currIdx) = secondHeDistance * vFactor;
        *(coeffTerms[0][2][0].first + currIdx) = secondVDistance;
        *(coeffTerms[0][2][1].first + currIdx) = secondVDistance * heFactor;
        *(coeffTerms[0][2][2].first + currIdx) = secondVDistance * vFactor;
        *(coeffTerms[1][1][0].first + currIdx) = firstHeDistance * secondHeDistance;
        *(coeffTerms[1][1][1].first + currIdx) = firstHeDistance * secondHeDistance * heFactor;
        *(coeffTerms[1][1][2].first + currIdx) = firstHeDistance * secondHeDistance * vFactor;
        *(coeffTerms[1][2][0].first + currIdx) = firstHeDistance * secondVDistance;
        *(coeffTerms[1][2][1].first + currIdx) = firstHeDistance * secondVDistance * heFactor;
        *(coeffTerms[1][2][2].first + currIdx) = firstHeDistance * secondVDistance * vFactor;
        *(coeffTerms[2][1][0].first + currIdx) = firstVDistance * secondHeDistance;
        *(coeffTerms[2][1][1].first + currIdx) = firstVDistance * secondHeDistance * heFactor;
        *(coeffTerms[2][1][2].first + currIdx) = firstVDistance * secondHeDistance * vFactor;
        *(coeffTerms[2][2][0].first + currIdx) = firstVDistance * secondVDistance;
        *(coeffTerms[2][2][1].first + currIdx) = firstVDistance * secondVDistance * heFactor;
        *(coeffTerms[2][2][2].first + currIdx) = firstVDistance * secondVDistance * vFactor;

        ++currIdx;
    });

    // Sum the terms to find deltas to the starting coefficient values,
    // and apply them.
    for(auto i = 0; i < 3; ++i)
    {
        for(auto j = 0; j < 3; ++j)
        {
            for(auto k = 0; k < 3; ++k)
            {
                auto currBeginIter = coeffTerms[i][j][k].first;
                auto currEndIter = coeffTerms[i][j][k].second;

                // Sort the values for the current coefficient.
                std::sort(currBeginIter, currEndIter);

                // Sum the values for the current coefficient.
                double currDelta = std::accumulate(currBeginIter, currEndIter, 0.0);

                // Update the corresponding coefficient in prodPair.
                prodPair.a[i][j][k] += currDelta;
            }
        }
    }

#else // defined(USE_SORTED_COEFF_SUMS)

    // Sum without sorting values.
    // Since we don't need to sort the elements used to compute
    // a given prodPair.a[i][j][k] value, we can compute them on the fly
    // and just update prodPair.a[i][j][k] as we go.
    // It is faster and more memory efficient, but has the potential
    // to incur larger numerical error in the computed coefficients than 
    // sorting and adding the values smallest to largest.
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&prodPair](const PendingProductionReactionInfo& currPRI) {

        // Use names corresponding to those in single version.
        int a = currPRI.numHe;
        int b = currPRI.numV;
        int c = currPRI.i;
        int d = currPRI.j;

        double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
                secondVDistance = 0.0;
        if (prodPair.first.getType() == ReactantType::PSISuper) {
            auto const& super = static_cast<PSICluster const&>(prodPair.first);
            firstHeDistance = super.getHeDistance(c);
            firstVDistance = super.getVDistance(d);
        }
        if (prodPair.second.getType() == ReactantType::PSISuper) {
            auto const& super = static_cast<PSICluster const&>(prodPair.second);
            secondHeDistance = super.getHeDistance(c);
            secondVDistance = super.getVDistance(d);
        }
        double heFactor = (double) (a - numHe) / dispersionHe;
        double vFactor = (double) (b - numV) / dispersionV;
        // First is A, second is B, in A + B -> this
        prodPair.a[0][0][0] += 1.0;
        prodPair.a[0][0][1] += heFactor;
        prodPair.a[0][0][2] += vFactor;
        prodPair.a[1][0][0] += firstHeDistance;
        prodPair.a[1][0][1] += firstHeDistance * heFactor;
        prodPair.a[1][0][2] += firstHeDistance * vFactor;
        prodPair.a[2][0][0] += firstVDistance;
        prodPair.a[2][0][1] += firstVDistance * heFactor;
        prodPair.a[2][0][2] += firstVDistance * vFactor;
        prodPair.a[0][1][0] += secondHeDistance;
        prodPair.a[0][1][1] += secondHeDistance * heFactor;
        prodPair.a[0][1][2] += secondHeDistance * vFactor;
        prodPair.a[0][2][0] += secondVDistance;
        prodPair.a[0][2][1] += secondVDistance * heFactor;
        prodPair.a[0][2][2] += secondVDistance * vFactor;
        prodPair.a[1][1][0] += firstHeDistance * secondHeDistance;
        prodPair.a[1][1][1] += firstHeDistance * secondHeDistance * heFactor;
        prodPair.a[1][1][2] += firstHeDistance * secondHeDistance * vFactor;
        prodPair.a[1][2][0] += firstHeDistance * secondVDistance;
        prodPair.a[1][2][1] += firstHeDistance * secondVDistance * heFactor;
        prodPair.a[1][2][2] += firstHeDistance * secondVDistance * vFactor;
        prodPair.a[2][1][0] += firstVDistance * secondHeDistance;
        prodPair.a[2][1][1] += firstVDistance * secondHeDistance * heFactor;
        prodPair.a[2][1][2] += firstVDistance * secondHeDistance * vFactor;
        prodPair.a[2][2][0] += firstVDistance * secondVDistance;
        prodPair.a[2][2][1] += firstVDistance * secondVDistance * heFactor;
        prodPair.a[2][2][2] += firstVDistance * secondVDistance * vFactor;
    });
#endif // defined(USE_SORTED_COEFF_SUMS)

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction, int a, int b) {

#if READY
    if (not partInProd_callsCounter) {
        partInProd_callsCounter = handlerRegistry->getEventCounter("PSISuper_partInProd_calls");
    }
    partInProd_callsCounter->increment();
#endif // READY

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>((reaction.first.getId() == id) ?
                                reaction.second :
                                reaction.first);

	// Check if we already know about the reaction.
    auto rkey = &otherCluster;
    auto it = effCombiningList.find(rkey);
	if (it == effCombiningList.end()) {

#if READY
        if (not partInProd_addsCounter) {
            partInProd_addsCounter = handlerRegistry->getEventCounter("PSISuper_partInProd_adds");
        }
        partInProd_addsCounter->increment();
#endif // READY

		// We did not already know about the reaction.
        // Note that we combine with the other cluster in this reaction.
        auto eret = effCombiningList.emplace(std::piecewise_construct,
                                std::forward_as_tuple(rkey),
                                std::forward_as_tuple(
                                    reaction,
                                    static_cast<PSICluster&>(otherCluster)));
        // Since we already checked and didn't know about the reaction then,
        // we had better have added it with our emplace call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effCombiningList.end());
    auto& combCluster = it->second;

	// Update the coefficients
	double heDistance = getHeDistance(a);
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vDistance = getVDistance(b);
	double vFactor = (double) (b - numV) / dispersionV;
	// This is A, itBis is B, in A + B -> C
	combCluster.a[0][0][0] += 1.0;
	combCluster.a[0][0][1] += heFactor;
	combCluster.a[0][0][2] += vFactor;
	combCluster.a[1][0][0] += heDistance;
	combCluster.a[1][0][1] += heDistance * heFactor;
	combCluster.a[1][0][2] += heDistance * vFactor;
	combCluster.a[2][0][0] += vDistance;
	combCluster.a[2][0][1] += vDistance * heFactor;
	combCluster.a[2][0][2] += vDistance * vFactor;

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction,
            const std::vector<PendingProductionReactionInfo>& prInfos) {

#if READY
    if (not partInProd_bcallsCounter) {
        partInProd_bcallsCounter = handlerRegistry->getEventCounter("PSISuper_partInProd_bcalls");
    }
    partInProd_bcallsCounter->increment();
#endif // READY

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>((reaction.first.getId() == id) ?
                                reaction.second :
                                reaction.first);

	// Check if we already know about the reaction.
    auto rkey = &otherCluster;
    auto it = effCombiningList.find(rkey);
	if (it == effCombiningList.end()) {

#if READY
        if (not partInProd_baddsCounter) {
            partInProd_baddsCounter = handlerRegistry->getEventCounter("PSISuper_partInProd_badds");
        }
        partInProd_baddsCounter->increment();
#endif // READY

		// We did not already know about the reaction.
        // Note that we combine with the other cluster in this reaction.
        auto eret = effCombiningList.emplace(std::piecewise_construct,
                                std::forward_as_tuple(rkey),
                                std::forward_as_tuple(
                                    reaction,
                                    static_cast<PSICluster&>(otherCluster)));
        // Since we already checked and didn't know about the reaction then,
        // we had better have added it with our emplace call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effCombiningList.end());
    auto& combCluster = it->second;

	// Update the coefficients
#if defined(USE_SORTED_COEFF_SUMS)

    // Update coefficients in a data-parallel way.
    // This requires more memory than the original approach since we 
    // store all of the values computed from the prInfos so that we can 
    // sort them before summing them for each (i,0,k) index in combCluster.a

    // Allocate space for all terms.
    // We do one large flat allocation, then drop "pointers" (iterators)
    // into it at appropriate offsets for each (i,0,k) triple.
    // This is more efficient than doing 9 memory allocations each
    // time we are called, and lets us still use the iterator-based C++
    // standard library functions like std::sort.
    std::vector<double> coeffTermsFlat(prInfos.size()*9);
    Array2D<std::pair<std::vector<double>::iterator, std::vector<double>::iterator>, 3, 3> coeffTerms;
    auto currIter = coeffTermsFlat.begin();
    for(auto i = 0; i < 3; ++i)
    {
        for(auto k = 0; k < 3; ++k)
        {
            auto nextIter = currIter + prInfos.size();
            coeffTerms[i][k] = std::make_pair(currIter, nextIter);
            currIter = nextIter;
        }
    }

    // Compute the values associated with each PRI.
    auto currIdx = 0;
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&combCluster,&coeffTerms,&currIdx](const PendingProductionReactionInfo& currPRI) {

            // Use names corresponding to the single-item version.
            int a = currPRI.i;
            int b = currPRI.j;

            double heDistance = getHeDistance(a);
            double heFactor = (double) (a - numHe) / dispersionHe;
            double vDistance = getVDistance(b);
            double vFactor = (double) (b - numV) / dispersionV;
            // This is A, itBis is B, in A + B -> C
            *(coeffTerms[0][0].first + currIdx) = 1.0;
            *(coeffTerms[0][1].first + currIdx) = heFactor;
            *(coeffTerms[0][2].first + currIdx) = vFactor;
            *(coeffTerms[1][0].first + currIdx) = heDistance;
            *(coeffTerms[1][1].first + currIdx) = heDistance * heFactor;
            *(coeffTerms[1][2].first + currIdx) = heDistance * vFactor;
            *(coeffTerms[2][0].first + currIdx) = vDistance;
            *(coeffTerms[2][1].first + currIdx) = vDistance * heFactor;
            *(coeffTerms[2][2].first + currIdx) = vDistance * vFactor;

            ++currIdx;
        });

    // Sum the terms to find deltas to the starting coefficient values,
    // and apply them.
    for(auto i = 0; i < 3; ++i)
    {
        for(auto k = 0; k < 3; ++k)
        {
            auto currBeginIter = coeffTerms[i][k].first;
            auto currEndIter = coeffTerms[i][k].second;

            // Sort the values for the current coefficient.
            std::sort(currBeginIter, currEndIter);

            // Sum the values for the current coefficient.
            double currDelta = std::accumulate(currBeginIter, currEndIter, 0.0);

            // Update the corresponding coefficient in prodPair.
            combCluster.a[i][0][k] += currDelta;
        }
    }
#else // defined(USE_SORTED_COEFF_SUMS)

    // Sum without sorting values.
    // Since we don't need to sort the elements used to compute
    // a given coefficient [i][j][k] value, we can compute them on the fly
    // and just update coefficient [i][j][k] as we go.
    // It is faster and more memory efficient, but has the potential
    // to incur larger numerical error in the computed coefficients than 
    // sorting and adding the values smallest to largest.
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&combCluster](const PendingProductionReactionInfo& currPRInfo) {

            // Use names corresponding to the single-item version.
            int a = currPRInfo.i;
            int b = currPRInfo.j;

            double heDistance = getHeDistance(a);
            double heFactor = (double) (a - numHe) / dispersionHe;
            double vDistance = getVDistance(b);
            double vFactor = (double) (b - numV) / dispersionV;
            // This is A, itBis is B, in A + B -> C
            combCluster.a[0][0][0] += 1.0;
            combCluster.a[0][0][1] += heFactor;
            combCluster.a[0][0][2] += vFactor;
            combCluster.a[1][0][0] += heDistance;
            combCluster.a[1][0][1] += heDistance * heFactor;
            combCluster.a[1][0][2] += heDistance * vFactor;
            combCluster.a[2][0][0] += vDistance;
            combCluster.a[2][0][1] += vDistance * heFactor;
            combCluster.a[2][0][2] += vDistance * vFactor;
        });
#endif // defined(USE_SORTED_COEFF_SUMS)

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction,
        int a, int b, int c, int d) {

#if READY
    if (not partInDiss_callsCounter) {
        partInDiss_callsCounter = handlerRegistry->getEventCounter("PSISuper_partInDiss_calls");
    }
    partInDiss_callsCounter->increment();
#endif // READY

	// Determine which is the other cluster.
	auto& emittedCluster = static_cast<PSICluster&>((reaction.first.getId() == id) ?
                            reaction.second :
                            reaction.first);

	// Check if we already know about the reaction.
    auto rkey = std::make_pair(&(reaction.dissociating), &emittedCluster);
    auto it = effDissociatingList.find(rkey);
	if (it == effDissociatingList.end()) {

		// We did not already know about it.
#if READY
        if (not partInDiss_addsCounter) {
            partInDiss_addsCounter = handlerRegistry->getEventCounter("PSISuper_partInDiss_adds");
        }
        partInDiss_addsCounter->increment();
#endif // READY

		// Add it to the network
        auto eret = effDissociatingList.emplace(std::piecewise_construct,
                        std::forward_as_tuple(rkey),
                        std::forward_as_tuple(
                            reaction,
				            static_cast<PSICluster&>(reaction.dissociating),
				            static_cast<PSICluster&>(emittedCluster)));
        // Since we already checked and didn't know about the reaction then,
        // we had better have added it with our emplace() call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effDissociatingList.end());
    auto& dissPair = it->second;

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSICluster const&>(reaction.dissociating);
		firstHeDistance = super.getHeDistance(a);
		firstVDistance = super.getVDistance(b);
	}
	double heFactor = (double) (c - numHe) / dispersionHe;
	double vFactor = (double) (d - numV) / dispersionV;

	// A is the dissociating cluster
	dissPair.a[0][0] += 1.0;
	dissPair.a[0][1] += heFactor;
	dissPair.a[0][2] += vFactor;
	dissPair.a[1][0] += firstHeDistance;
	dissPair.a[1][1] += firstHeDistance * heFactor;
	dissPair.a[1][2] += firstHeDistance * vFactor;
	dissPair.a[2][0] += firstVDistance;
	dissPair.a[2][1] += firstVDistance * heFactor;
	dissPair.a[2][2] += firstVDistance * vFactor;

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction,
        const std::vector<PendingProductionReactionInfo>& prInfos) {

#if READY
    if (not partInDiss_bcallsCounter) {
        partInDiss_bcallsCounter = handlerRegistry->getEventCounter("PSISuper_partInDiss_bcalls");
    }
    partInDiss_bcallsCounter->increment();
#endif // READY

	// Determine which is the other cluster.
	auto& emittedCluster = static_cast<PSICluster&>((reaction.first.getId() == id) ?
                            reaction.second :
                            reaction.first);

	// Check if we already know about the reaction.
    auto rkey = std::make_pair(&(reaction.dissociating), &emittedCluster);
    auto it = effDissociatingList.find(rkey);
	if (it == effDissociatingList.end()) {

		// We did not already know about it.
#if READY
        if (not partInDiss_baddsCounter) {
            partInDiss_baddsCounter = handlerRegistry->getEventCounter("PSISuper_partInDiss_badds");
        }
        partInDiss_baddsCounter->increment();
#endif // READY

		// Add it to the network
        auto eret = effDissociatingList.emplace(std::piecewise_construct,
                        std::forward_as_tuple(rkey),
                        std::forward_as_tuple(
                            reaction,
				            static_cast<PSICluster&>(reaction.dissociating),
				            static_cast<PSICluster&>(emittedCluster)));
        // Since we already checked and didn't know about the reaction then,
        // we had better have added it with our emplace() call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effDissociatingList.end());
    auto& dissPair = it->second;

	// Update the coefficients
#if defined(USE_SORTED_COEFF_SUMS)
    // Update coefficients in a data-parallel way.
    // This requires more memory than the original approach since we 
    // store all of the values computed from the prInfos so that we can 
    // sort them before summing them for each (i,j) value.

    // Allocate space for all terms.
    // We do one large flat allocation, then drop "pointers" (iterators)
    // into it at appropriate offsets for each (i,j) triple.
    // This is more efficient than doing 9 memory allocations each
    // time we are called, and lets us still use the iterator-based C++
    // standard library functions like std::sort.
    std::vector<double> coeffTermsFlat(prInfos.size()*9);
    Array2D<std::pair<std::vector<double>::iterator, std::vector<double>::iterator>, 3, 3> coeffTerms;
    auto currIter = coeffTermsFlat.begin();
    for(auto i = 0; i < 3; ++i)
    {
        for(auto j = 0; j < 3; ++j)
        {
            auto nextIter = currIter + prInfos.size();
            coeffTerms[i][j] = std::make_pair(currIter, nextIter);
            currIter = nextIter;
        }
    }

    auto currIdx = 0;
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&dissPair,&reaction,&coeffTerms,&currIdx](const PendingProductionReactionInfo& currPRI) {

            // Use names corresponding to the single-item version.
            int a = currPRI.numHe;
            int b = currPRI.numV;
            int c = currPRI.i;
            int d = currPRI.j;

            double firstHeDistance = 0.0, firstVDistance = 0.0;
            if (reaction.dissociating.getType() == ReactantType::PSISuper) {
                auto const& super = static_cast<PSICluster const&>(reaction.dissociating);
                firstHeDistance = super.getHeDistance(a);
                firstVDistance = super.getVDistance(b);
            }
            double heFactor = (double) (c - numHe) / dispersionHe;
            double vFactor = (double) (d - numV) / dispersionV;

            // A is the dissociating cluster
            *(coeffTerms[0][0].first + currIdx) = 1.0;
            *(coeffTerms[0][1].first + currIdx) = heFactor;
            *(coeffTerms[0][2].first + currIdx) = vFactor;
            *(coeffTerms[1][0].first + currIdx) = firstHeDistance;
            *(coeffTerms[1][1].first + currIdx) = firstHeDistance * heFactor;
            *(coeffTerms[1][2].first + currIdx) = firstHeDistance * vFactor;
            *(coeffTerms[2][0].first + currIdx) = firstVDistance;
            *(coeffTerms[2][1].first + currIdx) = firstVDistance * heFactor;
            *(coeffTerms[2][2].first + currIdx) = firstVDistance * vFactor;

            ++currIdx;
        });

    // Sum the terms to find deltas to the starting coefficient values,
    // and apply them.
    for(auto i = 0; i < 3; ++i)
    {
        for(auto j = 0; j < 3; ++j)
        {
            auto currBeginIter = coeffTerms[i][j].first;
            auto currEndIter = coeffTerms[i][j].second;

            // Sort the values for the current coefficient.
            std::sort(currBeginIter, currEndIter);

            // Sum the values for the current coefficient.
            double currDelta = std::accumulate(currBeginIter, currEndIter, 0.0);

            // Update the corresponding coefficient in prodPair.
            dissPair.a[i][j] += currDelta;
        }
    }
#else // defined(USE_SORTED_COEFF_SUMS)
    // Sum without sorting values.
    // Since we don't need to sort the elements used to compute
    // a given coefficient [i][j][k] value, we can compute them on the fly
    // and just update coefficient [i][j][k] as we go.
    // It is faster and more memory efficient, but has the potential
    // to incur larger numerical error in the computed coefficients than 
    // sorting and adding the values smallest to largest.
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&dissPair,&reaction](const PendingProductionReactionInfo& currPRI) {

            // Use names corresponding to the single-item version.
            int a = currPRI.numHe;
            int b = currPRI.numV;
            int c = currPRI.i;
            int d = currPRI.j;

            double firstHeDistance = 0.0, firstVDistance = 0.0;
            if (reaction.dissociating.getType() == ReactantType::PSISuper) {
                auto const& super = static_cast<PSICluster const&>(reaction.dissociating);
                firstHeDistance = super.getHeDistance(a);
                firstVDistance = super.getVDistance(b);
            }
            double heFactor = (double) (c - numHe) / dispersionHe;
            double vFactor = (double) (d - numV) / dispersionV;

            // A is the dissociating cluster
            dissPair.a[0][0] += 1.0;
            dissPair.a[0][1] += heFactor;
            dissPair.a[0][2] += vFactor;
            dissPair.a[1][0] += firstHeDistance;
            dissPair.a[1][1] += firstHeDistance * heFactor;
            dissPair.a[1][2] += firstHeDistance * vFactor;
            dissPair.a[2][0] += firstVDistance;
            dissPair.a[2][1] += firstVDistance * heFactor;
            dissPair.a[2][2] += firstVDistance * vFactor;
        });
#endif // defined(USE_SORTED_COEFF_SUMS)

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction,
        int a, int b, int c, int d) {

#if READY
    if (not emitFrom_callsCounter) {
        emitFrom_callsCounter = handlerRegistry->getEventCounter("PSICluster_emitFrom_calls");
    }
    emitFrom_callsCounter->increment();
#endif // READY

	// Check if we already know about the reaction.
    auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
    auto it = effEmissionList.find(rkey);
	if (it == effEmissionList.end()) {

		// We did not already know about it.
#if READY
        if (not emitFrom_addsCounter) {
            emitFrom_addsCounter = handlerRegistry->getEventCounter("PSICluster_emitFrom_adds");
        }
        emitFrom_addsCounter->increment();
#endif // READY

        // Note that we emit from the two rectants according to the given
        // reaction.
        auto eret = effEmissionList.emplace(std::piecewise_construct,
                        std::forward_as_tuple(rkey),
                        std::forward_as_tuple(
                            reaction,
                            static_cast<PSICluster&>(reaction.first),
                            static_cast<PSICluster&>(reaction.second)));
        // Since we already checked and didn't know about the reaction then,
        // we had better have added it with our emplace() call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effEmissionList.end());
    auto& dissPair = it->second;

	// Update the coefficients
	double heDistance = getHeDistance(a);
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vDistance = getVDistance(b);
	double vFactor = (double) (b - numV) / dispersionV;
	// A is the dissociating cluster
	dissPair.a[0][0] += 1.0;
	dissPair.a[0][1] += heFactor;
	dissPair.a[0][2] += vFactor;
	dissPair.a[1][0] += heDistance;
	dissPair.a[1][1] += heDistance * heFactor;
	dissPair.a[1][2] += heDistance * vFactor;
	dissPair.a[2][0] += vDistance;
	dissPair.a[2][1] += vDistance * heFactor;
	dissPair.a[2][2] += vDistance * vFactor;

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction,
        const std::vector<PendingProductionReactionInfo>& prInfos) {

#if READY
    if (not emitFrom_bcallsCounter) {
        emitFrom_bcallsCounter = handlerRegistry->getEventCounter("PSICluster_emitFrom_bcalls");
    }
    emitFrom_bcallsCounter->increment();
#endif // READY

	// Check if we already know about the reaction.
    auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
    auto it = effEmissionList.find(rkey);
	if (it == effEmissionList.end()) {

		// We did not already know about it.
#if READY
        if (not emitFrom_baddsCounter) {
            emitFrom_baddsCounter = handlerRegistry->getEventCounter("PSICluster_emitFrom_badds");
        }
        emitFrom_baddsCounter->increment();
#endif // READY

        // Note that we emit from the two rectants according to the given
        // reaction.
        auto eret = effEmissionList.emplace(std::piecewise_construct,
                        std::forward_as_tuple(rkey),
                        std::forward_as_tuple(
                            reaction,
                            static_cast<PSICluster&>(reaction.first),
                            static_cast<PSICluster&>(reaction.second)));
        // Since we already checked and didn't know about the reaction then,
        // we had better have added it with our emplace() call.
        assert(eret.second);
        it = eret.first;
	}
    assert(it != effEmissionList.end());
    auto& dissPair = it->second;

	// Update the coeeficients
#if defined(USE_SORTED_COEFF_SUMS)
    // Update coefficients in a data-parallel way.
    // This requires more memory than the original approach since we 
    // store all of the values computed from the prInfos so that we can 
    // sort them before summing them for each (i,j) value.

    // Allocate space for all terms.
    // We do one large flat allocation, then drop "pointers" (iterators)
    // into it at appropriate offsets for each (i,j) triple.
    // This is more efficient than doing 9 memory allocations each
    // time we are called, and lets us still use the iterator-based C++
    // standard library functions like std::sort.
    std::vector<double> coeffTermsFlat(prInfos.size()*9);
    Array2D<std::pair<std::vector<double>::iterator, std::vector<double>::iterator>, 3, 3> coeffTerms;
    auto currIter = coeffTermsFlat.begin();
    for(auto i = 0; i < 3; ++i)
    {
        for(auto j = 0; j < 3; ++j)
        {
            auto nextIter = currIter + prInfos.size();
            coeffTerms[i][j] = std::make_pair(currIter, nextIter);
            currIter = nextIter;
        }
    }

    auto currIdx = 0;
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&dissPair,&reaction,&coeffTerms,&currIdx](const PendingProductionReactionInfo& currPRI) {

            // Use same names as used in single version.
            int a = currPRI.numHe;
            int b = currPRI.numV;

            double heDistance = getHeDistance(a);
            double heFactor = (double) (a - numHe) / dispersionHe;
            double vDistance = getVDistance(b);
            double vFactor = (double) (b - numV) / dispersionV;
            // A is the dissociating cluster
            *(coeffTerms[0][0].first + currIdx) = 1.0;
            *(coeffTerms[0][1].first + currIdx) = heFactor;
            *(coeffTerms[0][2].first + currIdx) = vFactor;
            *(coeffTerms[1][0].first + currIdx) = heDistance;
            *(coeffTerms[1][1].first + currIdx) = heDistance * heFactor;
            *(coeffTerms[1][2].first + currIdx) = heDistance * vFactor;
            *(coeffTerms[2][0].first + currIdx) = vDistance;
            *(coeffTerms[2][1].first + currIdx) = vDistance * heFactor;
            *(coeffTerms[2][2].first + currIdx) = vDistance * vFactor;

            ++currIdx;
        });

    // Sum the terms to find deltas to the starting coefficient values,
    // and apply them.
    for(auto i = 0; i < 3; ++i)
    {
        for(auto j = 0; j < 3; ++j)
        {
            auto currBeginIter = coeffTerms[i][j].first;
            auto currEndIter = coeffTerms[i][j].second;

            // Sort the values for the current coefficient.
            std::sort(currBeginIter, currEndIter);

            // Sum the values for the current coefficient.
            double currDelta = std::accumulate(currBeginIter, currEndIter, 0.0);

            // Update the corresponding coefficient in prodPair.
            dissPair.a[i][j] += currDelta;
        }
    }
#else // defined(USE_SORTED_COEFF_SUMS)
    // Sum without sorting values.
    // Since we don't need to sort the elements used to compute
    // a given coefficient [i][j] value, we can compute them on the fly
    // and just update coefficient [i][j] as we go.
    // It is faster and more memory efficient, but has the potential
    // to incur larger numerical error in the computed coefficients than 
    // sorting and adding the values smallest to largest.
    std::for_each(prInfos.begin(), prInfos.end(),
        [this,&dissPair](const PendingProductionReactionInfo& currPRI) {

        // Use same names as used in single version.
        int a = currPRI.numHe;
        int b = currPRI.numV;

        double heDistance = getHeDistance(a);
        double heFactor = (double) (a - numHe) / dispersionHe;
        double vDistance = getVDistance(b);
        double vFactor = (double) (b - numV) / dispersionV;
        // A is the dissociating cluster
        dissPair.a[0][0] += 1.0;
        dissPair.a[0][1] += heFactor;
        dissPair.a[0][2] += vFactor;
        dissPair.a[1][0] += heDistance;
        dissPair.a[1][1] += heDistance * heFactor;
        dissPair.a[1][2] += heDistance * vFactor;
        dissPair.a[2][0] += vDistance;
        dissPair.a[2][1] += vDistance * heFactor;
        dissPair.a[2][2] += vDistance * vFactor;
    });
#endif // defined(USE_SORTED_COEFF_SUMS)

	return;
}

void PSISuperCluster::setHeVVector(std::vector<std::pair<int, int> > vec) {
	// Initialize the dispersion sum
	double nHeSquare = 0.0, nVSquare = 0.0;
	// Update the network map, compute the radius and dispersions
	for (auto it = vec.begin(); it != vec.end(); it++) {
		reactionRadius += xolotlCore::tungstenLatticeConstant
				* pow((3.0 * (double) ((*it).second)) / xolotlCore::pi,
						(1.0 / 3.0)) * 0.5 / (double) nTot;

		// Compute nSquare for the dispersion
		nHeSquare += (double) (*it).first * (*it).first;
		nVSquare += (double) (*it).second * (*it).second;
	}

	// Compute the dispersions
	if (sectionHeWidth == 1)
		dispersionHe = 1.0;
	else
		dispersionHe = 2.0 * (nHeSquare - (numHe * (double) nTot * numHe))
				/ ((double) (nTot * (sectionHeWidth - 1)));

	if (sectionVWidth == 1)
		dispersionV = 1.0;
	else
		dispersionV = 2.0 * (nVSquare - (numV * (double) nTot * numV))
				/ ((double) (nTot * (sectionVWidth - 1)));

	// Set the boundaries
    heBounds = IntegerRange<IReactant::SizeType>(
static_cast<IReactant::SizeType>((numHe - (double) sectionHeWidth / 2.0) + 1),
static_cast<IReactant::SizeType>((numHe - (double) sectionHeWidth / 2.0) + sectionHeWidth) + 1);
	vBounds = IntegerRange<IReactant::SizeType>(
static_cast<IReactant::SizeType>((numV - (double) sectionVWidth / 2.0) + 1),
static_cast<IReactant::SizeType>((numV - (double) sectionVWidth / 2.0) + sectionVWidth) + 1);

	return;
}

double PSISuperCluster::getTotalConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
    for (auto const& i : heBounds) {
        for (auto const& j : vBounds) {
			// Compute the distances
			heDistance = getHeDistance(i);
			vDistance = getVDistance(j);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance);
		}
	}

	return conc;
}

double PSISuperCluster::getTotalHeliumConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
    for (auto const& i : heBounds) {
        for (auto const& j : vBounds) {
			// Compute the distances
			heDistance = getHeDistance(i);
			vDistance = getVDistance(j);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) i;
		}
	}

	return conc;
}

double PSISuperCluster::getTotalVacancyConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
    for (auto const& i : heBounds) {
        for (auto const& j : vBounds) {
			// Compute the distances
			heDistance = getHeDistance(i);
			vDistance = getVDistance(j);

			// Add the concentration of each cluster in the group times its number of helium
			conc += getConcentration(heDistance, vDistance) * (double) j;
		}
	}

	return conc;
}

void PSISuperCluster::resetConnectivities() {
	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(heMomId);
	setDissociationConnectivity(heMomId);
	setReactionConnectivity(vMomId);
	setDissociationConnectivity(vMomId);

	// Visit all the reacting pairs
    std::for_each(effReactingList.begin(), effReactingList.end(),
        [this](ProductionPairMap::value_type const& currMapItem) {
            // The cluster is connecting to both clusters in the pair
            auto const& currPair = currMapItem.second;
            setReactionConnectivity(currPair.first.getId());
            setReactionConnectivity(currPair.first.getHeMomentumId());
            setReactionConnectivity(currPair.first.getVMomentumId());
            setReactionConnectivity(currPair.second.getId());
            setReactionConnectivity(currPair.second.getHeMomentumId());
            setReactionConnectivity(currPair.second.getVMomentumId());
        });

	// Visit all the combining pairs
    std::for_each(effCombiningList.begin(), effCombiningList.end(),
        [this](CombiningClusterMap::value_type const& currMapItem) {
            // The cluster is connecting to the combining cluster
            auto const& currComb = currMapItem.second;
            setReactionConnectivity(currComb.first.getId());
            setReactionConnectivity(currComb.first.getHeMomentumId());
            setReactionConnectivity(currComb.first.getVMomentumId());
        });

	// Loop over all the dissociating pairs
    std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
        [this](DissociationPairMap::value_type const& currMapItem) {
            // The cluster is connecting to the combining cluster
            auto const& currPair = currMapItem.second;
            setDissociationConnectivity(currPair.first.getId());
            setDissociationConnectivity(currPair.first.getHeMomentumId());
            setDissociationConnectivity(currPair.first.getVMomentumId());
        });

	// Don't loop on the effective emission pairs because
	// this cluster is not connected to them

	// Initialize the partial vector for the momentum
	int dof = network.getDOF();
	heMomentumPartials.resize(dof, 0.0);
	vMomentumPartials.resize(dof, 0.0);

	return;
}

double PSISuperCluster::getDissociationFlux() {
	// Initial declarations
	double flux = 0.0;

	// Sum over all the dissociating pairs
    // TODO consider using std::accumulate.  May also want to change side
    // effect of updating member variables heMomentumFlux and 
    // vMomentumFlux here.
    std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
        [this,&flux](DissociationPairMap::value_type const& currMapItem) {
            auto const& currPair = currMapItem.second;

            // Get the dissociating clusters
            auto const& dissociatingCluster = currPair.first;
            double l0A = dissociatingCluster.getConcentration(0.0, 0.0);
            double lHeA = dissociatingCluster.getHeMomentum();
            double lVA = dissociatingCluster.getVMomentum();
            // Update the flux
            auto value = currPair.kConstant / (double) nTot;
            flux += value * (currPair.a[0][0] * l0A + currPair.a[1][0] * lHeA + currPair.a[2][0] * lVA);
            // Compute the momentum fluxes
            heMomentumFlux += value
                    * (currPair.a[0][1] * l0A + currPair.a[1][1] * lHeA + currPair.a[2][1] * lVA);
            vMomentumFlux += value
                    * (currPair.a[0][2] * l0A + currPair.a[1][2] * lHeA + currPair.a[2][2] * lVA);
        });

	// Return the flux
	return flux;
}

double PSISuperCluster::getEmissionFlux() {
	// Initial declarations
	double flux = 0.0;

	// Loop over all the emission pairs
    // TODO consider using std::accumulate.  May also want to change side
    // effect of updating member variables heMomentumFlux and 
    // vMomentumFlux here.
    std::for_each(effEmissionList.begin(), effEmissionList.end(),
        [this,&flux](DissociationPairMap::value_type const& currMapItem) {
            auto const& currPair = currMapItem.second;

            // Update the flux
            auto value = currPair.kConstant / (double) nTot;
            flux += value * (currPair.a[0][0] * l0 + currPair.a[1][0] * l1He + currPair.a[2][0] * l1V);
            // Compute the momentum fluxes
            heMomentumFlux -= value
                    * (currPair.a[0][1] * l0 + currPair.a[1][1] * l1He + currPair.a[2][1] * l1V);
            vMomentumFlux -= value
                    * (currPair.a[0][2] * l0 + currPair.a[1][2] * l1He + currPair.a[2][2] * l1V);
        });

	return flux;
}

double PSISuperCluster::getProductionFlux() {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
    // TODO consider using std::accumulate.  May also want to change side
    // effect of updating member variables heMomentumFlux and 
    // vMomentumFlux here.
    std::for_each(effReactingList.begin(), effReactingList.end(),
        [this,&flux](ProductionPairMap::value_type const& currMapItem) {

        auto const& currPair = currMapItem.second;

		// Get the two reacting clusters
		auto const& firstReactant = currPair.first;
		auto const& secondReactant = currPair.second;
		double l0A = firstReactant.getConcentration(0.0, 0.0);
		double l0B = secondReactant.getConcentration(0.0, 0.0);
		double lHeA = firstReactant.getHeMomentum();
		double lHeB = secondReactant.getHeMomentum();
		double lVA = firstReactant.getVMomentum();
		double lVB = secondReactant.getVMomentum();
		// Update the flux
		auto value = currPair.kConstant / (double) nTot;
		flux += value
				* (currPair.a[0][0][0] * l0A * l0B + currPair.a[0][1][0] * l0A * lHeB
						+ currPair.a[0][2][0] * l0A * lVB + currPair.a[1][0][0] * lHeA * l0B
						+ currPair.a[1][1][0] * lHeA * lHeB + currPair.a[1][2][0] * lHeA * lVB
						+ currPair.a[2][0][0] * lVA * l0B + currPair.a[2][1][0] * lVA * lHeB
						+ currPair.a[2][2][0] * lVA * lVB);
		// Compute the momentum fluxes
		heMomentumFlux += value
				* (currPair.a[0][0][1] * l0A * l0B + currPair.a[0][1][1] * l0A * lHeB
						+ currPair.a[0][2][1] * l0A * lVB + currPair.a[1][0][1] * lHeA * l0B
						+ currPair.a[1][1][1] * lHeA * lHeB + currPair.a[1][2][1] * lHeA * lVB
						+ currPair.a[2][0][1] * lVA * l0B + currPair.a[2][1][1] * lVA * lHeB
						+ currPair.a[2][2][1] * lVA * lVB);
		vMomentumFlux += value
				* (currPair.a[0][0][2] * l0A * l0B + currPair.a[0][1][2] * l0A * lHeB
						+ currPair.a[0][2][2] * l0A * lVB + currPair.a[1][0][2] * lHeA * l0B
						+ currPair.a[1][1][2] * lHeA * lHeB + currPair.a[1][2][2] * lHeA * lVB
						+ currPair.a[2][0][2] * lVA * l0B + currPair.a[2][1][2] * lVA * lHeB
						+ currPair.a[2][2][2] * lVA * lVB);
	});

	// Return the production flux
	return flux;
}

double PSISuperCluster::getCombinationFlux() {
	// Local declarations
	double flux = 0.0;

	// Sum over all the combining clusters
    // TODO consider using std::accumulate.  May also want to change side
    // effect of updating member variables heMomentumFlux and 
    // vMomentumFlux here.
    std::for_each(effCombiningList.begin(), effCombiningList.end(),
        [this,&flux](CombiningClusterMap::value_type const& currMapItem) {
            // Get the combining cluster
            auto const& currComb = currMapItem.second;
            auto const& combiningCluster = currComb.first;
            double l0B = combiningCluster.getConcentration(0.0, 0.0);
            double lHeB = combiningCluster.getHeMomentum();
            double lVB = combiningCluster.getVMomentum();
            // Update the flux
            auto value = currComb.kConstant / (double) nTot;
            flux += value
                    * (currComb.a[0][0][0] * l0B * l0 + currComb.a[1][0][0] * l0B * l1He
                            + currComb.a[2][0][0] * l0B * l1V + currComb.a[0][1][0] * lHeB * l0
                            + currComb.a[1][1][0] * lHeB * l1He + currComb.a[2][1][0] * lHeB * l1V
                            + currComb.a[0][2][0] * lVB * l0 + currComb.a[1][2][0] * lVB * l1He
                            + currComb.a[2][2][0] * lVB * l1V);
            // Compute the momentum fluxes
            heMomentumFlux -= value
                    * (currComb.a[0][0][1] * l0B * l0 + currComb.a[1][0][1] * l0B * l1He
                            + currComb.a[2][0][1] * l0B * l1V + currComb.a[0][1][1] * lHeB * l0
                            + currComb.a[1][1][1] * lHeB * l1He + currComb.a[2][1][1] * lHeB * l1V
                            + currComb.a[0][2][1] * lVB * l0 + currComb.a[1][2][1] * lVB * l1He
                            + currComb.a[2][2][1] * lVB * l1V);
            vMomentumFlux -= value
                    * (currComb.a[0][0][2] * l0B * l0 + currComb.a[1][0][2] * l0B * l1He
                            + currComb.a[2][0][2] * l0B * l1V + currComb.a[0][1][2] * lHeB * l0
                            + currComb.a[1][1][2] * lHeB * l1He + currComb.a[2][1][2] * lHeB * l1V
                            + currComb.a[0][2][2] * lVB * l0 + currComb.a[1][2][2] * lVB * l1He
                            + currComb.a[2][2][2] * lVB * l1V);
        });

	return flux;
}

void PSISuperCluster::getPartialDerivatives(
		std::vector<double> & partials) const {
	// Reinitialize the momentum partial derivatives vector
	std::fill(heMomentumPartials.begin(), heMomentumPartials.end(), 0.0);
	std::fill(vMomentumPartials.begin(), vMomentumPartials.end(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void PSISuperCluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
    std::for_each(effReactingList.begin(), effReactingList.end(),
        [this,&partials](ProductionPairMap::value_type const& currMapItem) {

        auto const& currPair = currMapItem.second;

		// Get the two reacting clusters
		auto const& firstReactant = currPair.first;
		auto const& secondReactant = currPair.second;
		double l0A = firstReactant.getConcentration(0.0, 0.0);
		double l0B = secondReactant.getConcentration(0.0, 0.0);
		double lHeA = firstReactant.getHeMomentum();
		double lHeB = secondReactant.getHeMomentum();
		double lVA = firstReactant.getVMomentum();
		double lVB = secondReactant.getVMomentum();

		// Compute the contribution from the first part of the reacting pair
		auto value = currPair.kConstant / (double) nTot;
		auto index = firstReactant.getId() - 1;
		partials[index] += value
				* (currPair.a[0][0][0] * l0B + currPair.a[0][1][0] * lHeB + currPair.a[0][2][0] * lVB);
		heMomentumPartials[index] += value
				* (currPair.a[0][0][1] * l0B + currPair.a[0][1][1] * lHeB + currPair.a[0][2][1] * lVB);
		vMomentumPartials[index] += value
				* (currPair.a[0][0][2] * l0B + currPair.a[0][1][2] * lHeB + currPair.a[0][2][2] * lVB);
		index = firstReactant.getHeMomentumId() - 1;
		partials[index] += value
				* (currPair.a[1][0][0] * l0B + currPair.a[1][1][0] * lHeB + currPair.a[1][2][0] * lVB);
		heMomentumPartials[index] += value
				* (currPair.a[1][0][1] * l0B + currPair.a[1][1][1] * lHeB + currPair.a[1][2][1] * lVB);
		vMomentumPartials[index] += value
				* (currPair.a[1][0][2] * l0B + currPair.a[1][1][2] * lHeB + currPair.a[1][2][2] * lVB);
		index = firstReactant.getVMomentumId() - 1;
		partials[index] += value
				* (currPair.a[2][0][0] * l0B + currPair.a[2][1][0] * lHeB + currPair.a[2][2][0] * lVB);
		heMomentumPartials[index] += value
				* (currPair.a[2][0][1] * l0B + currPair.a[2][1][1] * lHeB + currPair.a[2][2][1] * lVB);
		vMomentumPartials[index] += value
				* (currPair.a[2][0][2] * l0B + currPair.a[2][1][2] * lHeB + currPair.a[2][2][2] * lVB);
		// Compute the contribution from the second part of the reacting pair
		index = secondReactant.getId() - 1;
		partials[index] += value
				* (currPair.a[0][0][0] * l0A + currPair.a[1][0][0] * lHeA + currPair.a[2][0][0] * lVA);
		heMomentumPartials[index] += value
				* (currPair.a[0][0][1] * l0A + currPair.a[1][0][1] * lHeA + currPair.a[2][0][1] * lVA);
		vMomentumPartials[index] += value
				* (currPair.a[0][0][2] * l0A + currPair.a[1][0][2] * lHeA + currPair.a[2][0][2] * lVA);
		index = secondReactant.getHeMomentumId() - 1;
		partials[index] += value
				* (currPair.a[0][1][0] * l0A + currPair.a[1][1][0] * lHeA + currPair.a[2][1][0] * lVA);
		heMomentumPartials[index] += value
				* (currPair.a[0][1][1] * l0A + currPair.a[1][1][1] * lHeA + currPair.a[2][1][1] * lVA);
		vMomentumPartials[index] += value
				* (currPair.a[0][1][2] * l0A + currPair.a[1][1][2] * lHeA + currPair.a[2][1][2] * lVA);
		index = secondReactant.getVMomentumId() - 1;
		partials[index] += value
				* (currPair.a[0][2][0] * l0A + currPair.a[1][2][0] * lHeA + currPair.a[2][2][0] * lVA);
		heMomentumPartials[index] += value
				* (currPair.a[0][2][1] * l0A + currPair.a[1][2][1] * lHeA + currPair.a[2][2][1] * lVA);
		vMomentumPartials[index] += value
				* (currPair.a[0][2][2] * l0A + currPair.a[1][2][2] * lHeA + currPair.a[2][2][2] * lVA);
	});

	return;
}

void PSISuperCluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Visit all the combining clusters
    std::for_each(effCombiningList.begin(), effCombiningList.end(),
        [this,&partials](CombiningClusterMap::value_type const& currMapItem) {
            // Get the combining clusters
            auto const& currComb = currMapItem.second;
            auto const& cluster = currComb.first;
            double l0B = cluster.getConcentration(0.0, 0.0);
            double lHeB = cluster.getHeMomentum();
            double lVB = cluster.getVMomentum();

            // Compute the contribution from the combining cluster
            auto value = currComb.kConstant / (double) nTot;
            auto index = cluster.getId() - 1;
            partials[index] -= value
                    * (currComb.a[0][0][0] * l0 + currComb.a[1][0][0] * l1He + currComb.a[2][0][0] * l1V);
            heMomentumPartials[index] -= value
                    * (currComb.a[0][0][1] * l0 + currComb.a[1][0][1] * l1He + currComb.a[2][0][1] * l1V);
            vMomentumPartials[index] -= value
                    * (currComb.a[0][0][2] * l0 + currComb.a[1][0][2] * l1He + currComb.a[2][0][2] * l1V);
            index = cluster.getHeMomentumId() - 1;
            partials[index] -= value
                    * (currComb.a[0][1][0] * l0 + currComb.a[1][1][0] * l1He + currComb.a[2][1][0] * l1V);
            heMomentumPartials[index] -= value
                    * (currComb.a[0][1][1] * l0 + currComb.a[1][1][1] * l1He + currComb.a[2][1][1] * l1V);
            vMomentumPartials[index] -= value
                    * (currComb.a[0][1][2] * l0 + currComb.a[1][1][2] * l1He + currComb.a[2][1][2] * l1V);
            index = cluster.getVMomentumId() - 1;
            partials[index] -= value
                    * (currComb.a[0][2][0] * l0 + currComb.a[1][2][0] * l1He + currComb.a[2][2][0] * l1V);
            heMomentumPartials[index] -= value
                    * (currComb.a[0][2][1] * l0 + currComb.a[1][2][1] * l1He + currComb.a[2][2][1] * l1V);
            vMomentumPartials[index] -= value
                    * (currComb.a[0][2][2] * l0 + currComb.a[1][2][2] * l1He + currComb.a[2][2][2] * l1V);
            // Compute the contribution from this cluster
            index = id - 1;
            partials[index] -= value
                    * (currComb.a[0][0][0] * l0B + currComb.a[0][1][0] * lHeB + currComb.a[0][2][0] * lVB);
            heMomentumPartials[index] -= value
                    * (currComb.a[0][0][1] * l0B + currComb.a[0][1][1] * lHeB + currComb.a[0][2][1] * lVB);
            vMomentumPartials[index] -= value
                    * (currComb.a[0][0][2] * l0B + currComb.a[0][1][2] * lHeB + currComb.a[0][2][2] * lVB);
            index = heMomId - 1;
            partials[index] -= value
                    * (currComb.a[1][0][0] * l0B + currComb.a[1][1][0] * lHeB + currComb.a[1][2][0] * lVB);
            heMomentumPartials[index] -= value
                    * (currComb.a[1][0][1] * l0B + currComb.a[1][1][1] * lHeB + currComb.a[1][2][1] * lVB);
            vMomentumPartials[index] -= value
                    * (currComb.a[1][0][2] * l0B + currComb.a[1][1][2] * lHeB + currComb.a[1][2][2] * lVB);
            index = vMomId - 1;
            partials[index] -= value
                    * (currComb.a[2][0][0] * l0B + currComb.a[2][1][0] * lHeB + currComb.a[2][2][0] * lVB);
            heMomentumPartials[index] -= value
                    * (currComb.a[2][0][1] * l0B + currComb.a[2][1][1] * lHeB + currComb.a[2][2][1] * lVB);
            vMomentumPartials[index] -= value
                    * (currComb.a[2][0][2] * l0B + currComb.a[2][1][2] * lHeB + currComb.a[2][2][2] * lVB);
        });

	return;
}

void PSISuperCluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Visit all the dissociating pairs
    std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
        [this,&partials](DissociationPairMap::value_type const& currMapItem) {
            auto& currPair = currMapItem.second;

            // Get the dissociating clusters
            auto const& cluster = currPair.first;
            // Compute the contribution from the dissociating cluster
            auto value = currPair.kConstant / (double) nTot;
            auto index = cluster.getId() - 1;
            partials[index] += value * (currPair.a[0][0]);
            heMomentumPartials[index] += value * (currPair.a[0][1]);
            vMomentumPartials[index] += value * (currPair.a[0][2]);
            index = cluster.getHeMomentumId() - 1;
            partials[index] += value * (currPair.a[1][0]);
            heMomentumPartials[index] += value * (currPair.a[1][1]);
            vMomentumPartials[index] += value * (currPair.a[1][2]);
            index = cluster.getVMomentumId() - 1;
            partials[index] += value * (currPair.a[2][0]);
            heMomentumPartials[index] += value * (currPair.a[2][1]);
            vMomentumPartials[index] += value * (currPair.a[2][2]);
        });

	return;
}

void PSISuperCluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Visit all the emission pairs
    std::for_each(effEmissionList.begin(), effEmissionList.end(),
        [this,&partials](DissociationPairMap::value_type const& currMapItem) {
            auto& currPair = currMapItem.second;

            // Compute the contribution from the dissociating cluster
            auto value = currPair.kConstant / (double) nTot;
            auto index = id - 1;
            partials[index] -= value * (currPair.a[0][0]);
            heMomentumPartials[index] -= value * (currPair.a[0][1]);
            vMomentumPartials[index] -= value * (currPair.a[0][2]);
            index = heMomId - 1;
            partials[index] -= value * (currPair.a[1][0]);
            heMomentumPartials[index] -= value * (currPair.a[1][1]);
            vMomentumPartials[index] -= value * (currPair.a[1][2]);
            index = vMomId - 1;
            partials[index] -= value * (currPair.a[2][0]);
            heMomentumPartials[index] -= value * (currPair.a[2][1]);
            vMomentumPartials[index] -= value * (currPair.a[2][2]);
        });

	return;
}

void PSISuperCluster::getHeMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = heMomentumPartials[i];
	}

	return;
}

void PSISuperCluster::getVMomentPartialDerivatives(
		std::vector<double> & partials) const {
	// Loop on the size of the vector
	for (int i = 0; i < partials.size(); i++) {
		// Set to the values that were already computed
		partials[i] = vMomentumPartials[i];
	}

	return;
}


void PSISuperCluster::dumpCoefficients(std::ostream& os, PSISuperCluster::ProductionCoefficientBase const& curr) const {

    os << "a[0-2][0-2][0-2]: ";
    for(auto& curr2D : curr.a) {
        for(auto& curr1D : curr2D) {
            std::copy(curr1D.begin(), curr1D.end(),
                std::ostream_iterator<double>(os, " "));
        }
    }
}

void PSISuperCluster::dumpCoefficients(std::ostream& os, PSISuperCluster::SuperClusterDissociationPair const& curr) const {

    os << "a[0-2][0-2]: ";
    for(auto& curr1D : curr.a) {
        std::copy(curr1D.begin(), curr1D.end(),
            std::ostream_iterator<double>(os, " "));
    }
}


void PSISuperCluster::outputCoefficientsTo(std::ostream& os) const {

    os << "id: " << id << '\n';
    os << "reacting: " << effReactingList.size() << '\n';
    std::for_each(effReactingList.begin(), effReactingList.end(),
    [this,&os](ProductionPairMap::value_type const& currMapItem) {
        auto const& currPair = currMapItem.second;
        os << "first: " << currPair.first.getId()
            << "; second: " << currPair.second.getId()
            << "; ";
        dumpCoefficients(os, currPair);
        os << '\n';
    });

    os << "combining: " << effCombiningList.size() << '\n';
    std::for_each(effCombiningList.begin(), effCombiningList.end(),
    [this,&os](CombiningClusterMap::value_type const& currMapItem) {
        auto const& currComb = currMapItem.second;
        os << "other: " << currComb.first.getId()
            << "; ";
        dumpCoefficients(os, currComb);
        os << '\n';
    });

    os << "dissociating: " << effDissociatingList.size() << '\n';
    std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
    [this,&os](DissociationPairMap::value_type const& currMapItem) {
        auto const& currPair = currMapItem.second;
        os << "first: " << currPair.first.getId()
            << "; second: " << currPair.second.getId()
            << "; ";
        dumpCoefficients(os, currPair);
        os << '\n';
    });

    os << "emitting: " << effEmissionList.size() << '\n';
    std::for_each(effEmissionList.begin(), effEmissionList.end(),
    [this,&os](DissociationPairMap::value_type const& currMapItem) {
        auto const& currPair = currMapItem.second;
        os << "first: " << currPair.first.getId()
            << "; second: " << currPair.second.getId()
            << "; ";
        dumpCoefficients(os, currPair);
        os << '\n';
    });
}


