// Includes
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

double PSISuperCluster::getPade(double distHe, double distV) const {
	// Check if we are using the approximation
	if (padeCoef.size() == 0)
		return 0.0;

	// Initial declaration
	double numerator = 0.0, denominator = 0.0;

	numerator = padeCoef[0] + padeCoef[1] * distHe + padeCoef[2] * distV
			+ padeCoef[3] * distHe * distHe + padeCoef[4] * distV * distV
			+ padeCoef[5] * distHe * distV
			+ padeCoef[6] * distHe * distHe * distHe
			+ padeCoef[7] * distV * distV * distV
			+ padeCoef[8] * distHe * distHe * distV
			+ padeCoef[9] * distHe * distV * distV;

	denominator = 1.0 + padeCoef[10] * distHe + padeCoef[11] * distV
			+ padeCoef[12] * distHe * distHe + padeCoef[13] * distV * distV
			+ padeCoef[14] * distHe * distV;

//	std::cout << name << " "
//			<< distHe << " " << distV << " "
//			<< numerator / denominator << std::endl;

	return numerator / denominator;
}

PSISuperCluster::PSISuperCluster(double _numHe, double _numV, int _nTot,
		int heWidth, int vWidth, IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(_numHe, _numV)), numHe(_numHe), numV(
				_numV), nTot(_nTot), l0(0.0), l1He(0.0), l1V(0.0), dispersionHe(
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

void PSISuperCluster::resultFrom(ProductionReaction& reaction, int a, int b,
		int c, int d) {

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
				std::forward_as_tuple(reaction,
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
			secondVDistance = 0.0, firstPade = 0.0, secondPade = 0.0;
	if (prodPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(prodPair.first);
		firstHeDistance = super.getHeDistance(c);
		firstVDistance = super.getVDistance(d);
		firstPade = super.getPade(firstHeDistance, firstVDistance);
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(prodPair.second);
		secondHeDistance = super.getHeDistance(c);
		secondVDistance = super.getVDistance(d);
		secondPade = super.getPade(secondHeDistance, secondVDistance);
	}
	double heFactor = (double) (a - numHe) / dispersionHe;
	double vFactor = (double) (b - numV) / dispersionV;
	// First is A, second is B, in A + B -> this
	prodPair.a000 += 1.0;
	prodPair.a001 += heFactor;
	prodPair.a002 += vFactor;
	prodPair.a100 += firstHeDistance;
	prodPair.a101 += firstHeDistance * heFactor;
	prodPair.a102 += firstHeDistance * vFactor;
	prodPair.a200 += firstVDistance;
	prodPair.a201 += firstVDistance * heFactor;
	prodPair.a202 += firstVDistance * vFactor;
	prodPair.a010 += secondHeDistance;
	prodPair.a011 += secondHeDistance * heFactor;
	prodPair.a012 += secondHeDistance * vFactor;
	prodPair.a020 += secondVDistance;
	prodPair.a021 += secondVDistance * heFactor;
	prodPair.a022 += secondVDistance * vFactor;
	prodPair.a110 += firstHeDistance * secondHeDistance;
	prodPair.a111 += firstHeDistance * secondHeDistance * heFactor;
	prodPair.a112 += firstHeDistance * secondHeDistance * vFactor;
	prodPair.a120 += firstHeDistance * secondVDistance;
	prodPair.a121 += firstHeDistance * secondVDistance * heFactor;
	prodPair.a122 += firstHeDistance * secondVDistance * vFactor;
	prodPair.a210 += firstVDistance * secondHeDistance;
	prodPair.a211 += firstVDistance * secondHeDistance * heFactor;
	prodPair.a212 += firstVDistance * secondHeDistance * vFactor;
	prodPair.a220 += firstVDistance * secondVDistance;
	prodPair.a221 += firstVDistance * secondVDistance * heFactor;
	prodPair.a222 += firstVDistance * secondVDistance * vFactor;
	prodPair.b0 += firstPade * secondPade;
	prodPair.b1 += secondPade;
	prodPair.b2 += firstPade;
	prodPair.b3 += secondPade * firstHeDistance;
	prodPair.b4 += secondPade * firstVDistance;
	prodPair.b5 += firstPade * secondHeDistance;
	prodPair.b6 += firstPade * secondVDistance;
	prodPair.c0 += firstPade * secondPade * heFactor;
	prodPair.c1 += secondPade * heFactor;
	prodPair.c2 += firstPade * heFactor;
	prodPair.c3 += secondPade * firstHeDistance * heFactor;
	prodPair.c4 += secondPade * firstVDistance * heFactor;
	prodPair.c5 += firstPade * secondHeDistance * heFactor;
	prodPair.c6 += firstPade * secondVDistance * heFactor;
	prodPair.d0 += firstPade * secondPade * vFactor;
	prodPair.d1 += secondPade * vFactor;
	prodPair.d2 += firstPade * vFactor;
	prodPair.d3 += secondPade * firstHeDistance * vFactor;
	prodPair.d4 += secondPade * firstVDistance * vFactor;
	prodPair.d5 += firstPade * secondHeDistance * vFactor;
	prodPair.d6 += firstPade * secondVDistance * vFactor;

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
				std::forward_as_tuple(reaction,
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
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&prodPair](const PendingProductionReactionInfo& currPRI) {

				// Use names corresponding to those in single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				double firstHeDistance = 0.0, firstVDistance = 0.0, secondHeDistance = 0.0,
				secondVDistance = 0.0, firstPade = 0.0, secondPade = 0.0;
				if (prodPair.first.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSISuperCluster const&>(prodPair.first);
					firstHeDistance = super.getHeDistance(c);
					firstVDistance = super.getVDistance(d);
					firstPade = super.getPade(firstHeDistance, firstVDistance);
				}
				if (prodPair.second.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSISuperCluster const&>(prodPair.second);
					secondHeDistance = super.getHeDistance(c);
					secondVDistance = super.getVDistance(d);
					secondPade = super.getPade(secondHeDistance, secondVDistance);
				}

				double heFactor = (double) (a - numHe) / dispersionHe;
				double vFactor = (double) (b - numV) / dispersionV;
				// First is A, second is B, in A + B -> this
				prodPair.a000 += 1.0;
				prodPair.a001 += heFactor;
				prodPair.a002 += vFactor;
				prodPair.a100 += firstHeDistance;
				prodPair.a101 += firstHeDistance * heFactor;
				prodPair.a102 += firstHeDistance * vFactor;
				prodPair.a200 += firstVDistance;
				prodPair.a201 += firstVDistance * heFactor;
				prodPair.a202 += firstVDistance * vFactor;
				prodPair.a010 += secondHeDistance;
				prodPair.a011 += secondHeDistance * heFactor;
				prodPair.a012 += secondHeDistance * vFactor;
				prodPair.a020 += secondVDistance;
				prodPair.a021 += secondVDistance * heFactor;
				prodPair.a022 += secondVDistance * vFactor;
				prodPair.a110 += firstHeDistance * secondHeDistance;
				prodPair.a111 += firstHeDistance * secondHeDistance * heFactor;
				prodPair.a112 += firstHeDistance * secondHeDistance * vFactor;
				prodPair.a120 += firstHeDistance * secondVDistance;
				prodPair.a121 += firstHeDistance * secondVDistance * heFactor;
				prodPair.a122 += firstHeDistance * secondVDistance * vFactor;
				prodPair.a210 += firstVDistance * secondHeDistance;
				prodPair.a211 += firstVDistance * secondHeDistance * heFactor;
				prodPair.a212 += firstVDistance * secondHeDistance * vFactor;
				prodPair.a220 += firstVDistance * secondVDistance;
				prodPair.a221 += firstVDistance * secondVDistance * heFactor;
				prodPair.a222 += firstVDistance * secondVDistance * vFactor;
				prodPair.b0 += firstPade * secondPade;
				prodPair.b1 += secondPade;
				prodPair.b2 += firstPade;
				prodPair.b3 += secondPade * firstHeDistance;
				prodPair.b4 += secondPade * firstVDistance;
				prodPair.b5 += firstPade * secondHeDistance;
				prodPair.b6 += firstPade * secondVDistance;
				prodPair.c0 += firstPade * secondPade * heFactor;
				prodPair.c1 += secondPade * heFactor;
				prodPair.c2 += firstPade * heFactor;
				prodPair.c3 += secondPade * firstHeDistance * heFactor;
				prodPair.c4 += secondPade * firstVDistance * heFactor;
				prodPair.c5 += firstPade * secondHeDistance * heFactor;
				prodPair.c6 += firstPade * secondVDistance * heFactor;
				prodPair.d0 += firstPade * secondPade * vFactor;
				prodPair.d1 += secondPade * vFactor;
				prodPair.d2 += firstPade * vFactor;
				prodPair.d3 += secondPade * firstHeDistance * vFactor;
				prodPair.d4 += secondPade * firstVDistance * vFactor;
				prodPair.d5 += firstPade * secondHeDistance * vFactor;
				prodPair.d6 += firstPade * secondVDistance * vFactor;
			});

//	if (fabs(prodPair.b2) > 1.0e-6 || fabs(prodPair.b1) > 1.0e-6)
//		std::cout << name << " " << prodPair.first.getName() << " + "
//				<< prodPair.second.getName() << ": " << prodPair.b2 << " "
//				<< prodPair.b1 << std::endl;
	return;
}

void PSISuperCluster::resultFrom(ProductionReaction& reaction,
		IReactant& product) {

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto it = effReactingList.find(rkey);
	if (it == effReactingList.end()) {

		// We did not already know about this reaction.
		// Add info about production to our list.
		auto eret = effReactingList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
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

	// Values for grouping parameters
	int productLoHe = *(heBounds.begin()), productHiHe = *(heBounds.end()) - 1,
			productLoV = *(vBounds.begin()), productHiV = *(vBounds.end()) - 1,
			loHe = 0, hiHe = 0, loV = 0, hiV = 0, singleHeSize = 0,
			singleVSize = 0;

	if (prodPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(prodPair.first);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = prodPair.second.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}
	if (prodPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(prodPair.second);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = prodPair.first.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}

	int heWidth = std::min(productHiHe, hiHe + singleHeSize)
			- std::max(productLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(productHiV, hiV + singleVSize)
			- std::max(productLoV, loV + singleVSize) + 1;

	prodPair.a000 = heWidth * vWidth;
	prodPair.a001 = ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(productLoHe, loHe + singleHeSize),
					std::min(productHiHe, hiHe + singleHeSize), numHe);
	prodPair.a002 = ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(productLoV, loV + singleVSize),
					std::min(productHiV, hiV + singleVSize), numV);

	double a10 = ((double) (2 * vWidth) / (double) (hiHe - loHe))
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe),
					(double) (loHe + hiHe) / 2.0);
	prodPair.a010 = prodPair.second.isMixed() * a10;
	prodPair.a100 = prodPair.first.isMixed() * a10;

	double a20 = ((double) (2 * heWidth) / (double) (hiV - loV))
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV),
					(double) (loV + hiV) / 2.0);
	prodPair.a020 = prodPair.second.isMixed() * a20;
	prodPair.a200 = prodPair.first.isMixed() * a20;

	double a11 = ((double) (2 * vWidth)
			/ ((double) (hiHe - loHe) * dispersionHe))
			* secondOrderOffsetSum(std::max(productLoHe, loHe + singleHeSize),
					std::min(productHiHe, hiHe + singleHeSize), numHe,
					(double) (loHe + hiHe) / 2.0, -singleHeSize);
	prodPair.a011 = prodPair.second.isMixed() * a11;
	prodPair.a101 = prodPair.first.isMixed() * a11;

	double a12 = ((double) (2 * vWidth) / (double) (hiV - loV))
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe),
					(double) (loHe + hiHe) / 2.0)
			* ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(productLoV, loV + singleVSize),
					std::min(productHiV, hiV + singleVSize), numV);
	prodPair.a012 = prodPair.second.isMixed() * a12;
	prodPair.a102 = prodPair.first.isMixed() * a12;

	double a21 = ((double) (2 * heWidth) / (double) (hiHe - loHe))
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV),
					(double) (loV + hiV) / 2.0)
			* ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(productLoHe, loHe + singleHeSize),
					std::min(productHiHe, hiHe + singleHeSize), numHe);
	prodPair.a021 = prodPair.second.isMixed() * a21;
	prodPair.a201 = prodPair.first.isMixed() * a21;

	double a22 = ((double) (2 * heWidth) / ((double) (hiV - loV) * dispersionV))
			* secondOrderOffsetSum(std::max(productLoV, loV + singleVSize),
					std::min(productHiV, hiV + singleVSize), numV,
					(double) (loV + hiV) / 2.0, -singleVSize);
	prodPair.a022 = prodPair.second.isMixed() * a22;
	prodPair.a202 = prodPair.first.isMixed() * a22;

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction, int a,
		int b) {

#if READY
	if (not partInProd_callsCounter) {
		partInProd_callsCounter = handlerRegistry->getEventCounter("PSISuper_partInProd_calls");
	}
	partInProd_callsCounter->increment();
#endif // READY

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

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
				std::forward_as_tuple(reaction,
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
	double pade = getPade(heDistance, vDistance);
	// This is A, itBis is B, in A + B -> C
	combCluster.a000 += 1.0;
	combCluster.a001 += heFactor;
	combCluster.a002 += vFactor;
	combCluster.a100 += heDistance;
	combCluster.a101 += heDistance * heFactor;
	combCluster.a102 += heDistance * vFactor;
	combCluster.a200 += vDistance;
	combCluster.a201 += vDistance * heFactor;
	combCluster.a202 += vDistance * vFactor;
	combCluster.b2 += pade;
	combCluster.c2 += pade * heFactor;
	combCluster.d2 += pade * vFactor;

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction,
		const std::vector<PendingProductionReactionInfo>& pendingPRInfos) {

#if READY
	if (not partInProd_bcallsCounter) {
		partInProd_bcallsCounter = handlerRegistry->getEventCounter("PSISuper_partInProd_bcalls");
	}
	partInProd_bcallsCounter->increment();
#endif // READY

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

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
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(otherCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effCombiningList.end());
	auto& combCluster = it->second;

	// Update the coefficients
	std::for_each(pendingPRInfos.begin(), pendingPRInfos.end(),
			[this,&combCluster,&otherCluster](const PendingProductionReactionInfo& currPRInfo) {

				// Use names corresponding to the single-item version.
				int a = currPRInfo.i;
				int b = currPRInfo.j;

				double heDistance = getHeDistance(a);
				double heFactor = (double) (a - numHe) / dispersionHe;
				double vDistance = getVDistance(b);
				double vFactor = (double) (b - numV) / dispersionV;
				double pade = getPade(heDistance, vDistance);
				// This is A, itBis is B, in A + B -> C
				combCluster.a000 += 1.0;
				combCluster.a001 += heFactor;
				combCluster.a002 += vFactor;
				combCluster.a100 += heDistance;
				combCluster.a101 += heDistance * heFactor;
				combCluster.a102 += heDistance * vFactor;
				combCluster.a200 += vDistance;
				combCluster.a201 += vDistance * heFactor;
				combCluster.a202 += vDistance * vFactor;
				combCluster.b2 += pade;
				combCluster.c2 += pade * heFactor;
				combCluster.d2 += pade * vFactor;
			});

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction,
		IReactant& prod) {

	setReactionConnectivity(id);
	// Look for the other cluster
	auto& otherCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if we already know about the reaction.
	auto rkey = &otherCluster;
	auto it = effCombiningList.find(rkey);
	if (it == effCombiningList.end()) {

		// We did not already know about the reaction.
		// Note that we combine with the other cluster in this reaction.
		auto eret = effCombiningList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(otherCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effCombiningList.end());
	auto& combCluster = it->second;

	// Values for grouping parameters
	int loHe = *(heBounds.begin()), hiHe = *(heBounds.end()) - 1, loV =
			*(vBounds.begin()), hiV = *(vBounds.end()) - 1, productLoHe = 0,
			productHiHe = 0, productLoV = 0, productHiV = 0;

	auto singleComp = otherCluster.getComposition();
	int singleHeSize = singleComp[toCompIdx(Species::He)];
	int singleVSize = singleComp[toCompIdx(Species::V)]
			- singleComp[toCompIdx(Species::I)];        // can be < 0

	if (prod.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(prod);
		auto const& heBounds = super.getHeBounds();
		productLoHe = *(heBounds.begin());
		productHiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		productLoV = *(vBounds.begin());
		productHiV = *(vBounds.end()) - 1;
	} else {
		auto productComp = prod.getComposition();
		productLoHe = productComp[toCompIdx(Species::He)];
		productHiHe = productComp[toCompIdx(Species::He)];
		productLoV = productComp[toCompIdx(Species::V)];
		productHiV = productComp[toCompIdx(Species::V)];
	}

	int heWidth = std::min(productHiHe, hiHe + singleHeSize)
			- std::max(productLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(productHiV, hiV + singleVSize)
			- std::max(productLoV, loV + singleVSize) + 1;

	combCluster.a000 += heWidth * vWidth;

	combCluster.a001 += ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe), numHe);

	combCluster.a002 += ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV), numV);

	combCluster.a100 += ((double) (2 * vWidth) / (double) (sectionHeWidth - 1))
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe), numHe);

	combCluster.a200 += ((double) (2 * heWidth) / (double) (sectionVWidth - 1))
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV), numV);

	combCluster.a101 += ((double) (2 * vWidth)
			/ ((double) (sectionHeWidth - 1) * dispersionHe))
			* secondOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe), numHe);

	combCluster.a102 += ((double) (2 * vWidth) / (double) (sectionHeWidth - 1))
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe), numHe)
			* ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV), numV);

	combCluster.a201 += ((double) (2 * heWidth) / (double) (sectionVWidth - 1))
			* firstOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV), numV)
			* ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(productLoHe - singleHeSize, loHe),
					std::min(productHiHe - singleHeSize, hiHe), numHe);

	combCluster.a202 += ((double) (2 * heWidth)
			/ ((double) (sectionVWidth - 1) * dispersionV))
			* secondOrderSum(std::max(productLoV - singleVSize, loV),
					std::min(productHiV - singleVSize, hiV), numV);

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction, int a,
		int b, int c, int d) {

#if READY
	if (not partInDiss_callsCounter) {
		partInDiss_callsCounter = handlerRegistry->getEventCounter("PSISuper_partInDiss_calls");
	}
	partInDiss_callsCounter->increment();
#endif // READY

	// Determine which is the other cluster.
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

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
				std::forward_as_tuple(reaction,
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
	double firstHeDistance = 0.0, firstVDistance = 0.0, firstPade = 0.0;
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super =
				static_cast<PSISuperCluster const&>(reaction.dissociating);
		firstHeDistance = super.getHeDistance(a);
		firstVDistance = super.getVDistance(b);
		firstPade = super.getPade(firstHeDistance, firstVDistance);
	}
	double heFactor = (double) (c - numHe) / dispersionHe;
	double vFactor = (double) (d - numV) / dispersionV;

	// A is the dissociating cluster
	dissPair.a00 += 1.0;
	dissPair.a01 += heFactor;
	dissPair.a02 += vFactor;
	dissPair.a10 += firstHeDistance;
	dissPair.a11 += firstHeDistance * heFactor;
	dissPair.a12 += firstHeDistance * vFactor;
	dissPair.a20 += firstVDistance;
	dissPair.a21 += firstVDistance * heFactor;
	dissPair.a22 += firstVDistance * vFactor;
	dissPair.b0 += firstPade;
	dissPair.c0 += firstPade * heFactor;
	dissPair.d0 += firstPade * vFactor;

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
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

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
				std::forward_as_tuple(reaction,
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
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&dissPair,&reaction](const PendingProductionReactionInfo& currPRI) {

				// Use names corresponding to the single-item version.
				int a = currPRI.numHe;
				int b = currPRI.numV;
				int c = currPRI.i;
				int d = currPRI.j;

				double firstHeDistance = 0.0, firstVDistance = 0.0, firstPade = 0.0;
				if (reaction.dissociating.getType() == ReactantType::PSISuper) {
					auto const& super = static_cast<PSISuperCluster const&>(reaction.dissociating);
					firstHeDistance = super.getHeDistance(a);
					firstVDistance = super.getVDistance(b);
					firstPade = super.getPade(firstHeDistance, firstVDistance);
				}
				double heFactor = (double) (c - numHe) / dispersionHe;
				double vFactor = (double) (d - numV) / dispersionV;

				// A is the dissociating cluster
				dissPair.a00 += 1.0;
				dissPair.a01 += heFactor;
				dissPair.a02 += vFactor;
				dissPair.a10 += firstHeDistance;
				dissPair.a11 += firstHeDistance * heFactor;
				dissPair.a12 += firstHeDistance * vFactor;
				dissPair.a20 += firstVDistance;
				dissPair.a21 += firstVDistance * heFactor;
				dissPair.a22 += firstVDistance * vFactor;
				dissPair.b0 += firstPade;
				dissPair.c0 += firstPade * heFactor;
				dissPair.d0 += firstPade * vFactor;
			});

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction,
		IReactant& disso) {

	// Determine which is the other cluster.
	auto& emittedCluster = static_cast<PSICluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.dissociating), &emittedCluster);
	auto it = effDissociatingList.find(rkey);
	if (it == effDissociatingList.end()) {

		// We did not already know about it.

		// Add it to the network
		auto eret = effDissociatingList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.dissociating),
						static_cast<PSICluster&>(emittedCluster)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effDissociatingList.end());
	auto& dissPair = it->second;

	// Values for grouping parameters
	int dissoLoHe = 0, dissoHiHe = 0, dissoLoV = 0, dissoHiV = 0, loHe =
			*(heBounds.begin()), hiHe = *(heBounds.end()) - 1, loV =
			*(vBounds.begin()), hiV = *(vBounds.end()) - 1;

	if (disso.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(disso);
		auto const& heBounds = super.getHeBounds();
		dissoLoHe = *(heBounds.begin());
		dissoHiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		dissoLoV = *(vBounds.begin());
		dissoHiV = *(vBounds.end()) - 1;
	} else {
		auto dissoComp = disso.getComposition();
		dissoLoHe = dissoComp[toCompIdx(Species::He)];
		dissoHiHe = dissoComp[toCompIdx(Species::He)];
		dissoLoV = dissoComp[toCompIdx(Species::V)];
		dissoHiV = dissoComp[toCompIdx(Species::V)];
	}
	auto singleComp = emittedCluster.getComposition();
	int singleHeSize = singleComp[toCompIdx(Species::He)];
	int singleVSize = singleComp[toCompIdx(Species::V)]
			- singleComp[toCompIdx(Species::I)]; // can be < 0

	int heWidth = std::min(dissoHiHe, hiHe + singleHeSize)
			- std::max(dissoLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(dissoHiV, hiV + singleVSize)
			- std::max(dissoLoV, loV + singleVSize) + 1;

	dissPair.a00 = heWidth * vWidth;

	dissPair.a01 = ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(dissoLoHe - singleHeSize, loHe),
					std::min(dissoHiHe - singleHeSize, hiHe), numHe);

	dissPair.a02 = ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(dissoLoV - singleVSize, loV),
					std::min(dissoHiV - singleVSize, hiV), numV);

	dissPair.a10 = ((double) (2 * vWidth) / (double) (dissoHiHe - dissoLoHe))
			* firstOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize),
					(dissoLoHe + dissoHiHe) / 2.0);

	dissPair.a20 = ((double) (2 * heWidth) / (double) (dissoHiV - dissoLoV))
			* firstOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize),
					(dissoLoV + dissoHiV) / 2.0);

	dissPair.a11 = ((double) (2 * vWidth)
			/ ((double) (dissoHiHe - dissoLoHe) * dispersionHe))
			* secondOrderOffsetSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize),
					(double) (dissoLoHe + dissoHiHe) / 2.0, numHe,
					-singleHeSize);

	dissPair.a12 = ((double) (2 * vWidth) / (double) (dissoHiV - dissoLoV))
			* firstOrderSum(std::max(dissoLoHe - singleHeSize, loHe),
					std::min(dissoHiHe - singleHeSize, hiHe), numHe)
			* ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(dissoLoV - singleVSize, loV),
					std::min(dissoHiV - singleVSize, hiV), numV);

	dissPair.a21 = ((double) (2 * heWidth) / (double) (dissoHiHe - dissoLoHe))
			* firstOrderSum(std::max(dissoLoV - singleVSize, loV),
					std::min(dissoHiV - singleVSize, hiV), numV)
			* ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(dissoLoHe - singleHeSize, loHe),
					std::min(dissoHiHe - singleHeSize, hiHe), numHe);

	dissPair.a22 = ((double) (2 * heWidth)
			/ ((double) (dissoHiV - dissoLoV) * dispersionV))
			* secondOrderOffsetSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize),
					(double) (dissoLoV + dissoHiV) / 2.0, numV, -singleVSize);

	if (disso.getType() != ReactantType::PSISuper) {
		dissPair.a10 = 0.0, dissPair.a20 = 0.0, dissPair.a11 = 0.0, dissPair.a12 =
				0.0, dissPair.a21 = 0.0, dissPair.a22 = 0.0;
	}

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction, int a, int b,
		int c, int d) {

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
				std::forward_as_tuple(reaction,
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
	double pade = getPade(heDistance, vDistance);
	// A is the dissociating cluster
	dissPair.a00 += 1.0;
	dissPair.a01 += heFactor;
	dissPair.a02 += vFactor;
	dissPair.a10 += heDistance;
	dissPair.a11 += heDistance * heFactor;
	dissPair.a12 += heDistance * vFactor;
	dissPair.a20 += vDistance;
	dissPair.a21 += vDistance * heFactor;
	dissPair.a22 += vDistance * vFactor;
	dissPair.b0 += pade;
	dissPair.c0 += pade * heFactor;
	dissPair.d0 += pade * vFactor;

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
				std::forward_as_tuple(reaction,
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
	std::for_each(prInfos.begin(), prInfos.end(),
			[this,&dissPair](const PendingProductionReactionInfo& currPRI) {

				// Use same names as used in single version.
				int a = currPRI.numHe;
				int b = currPRI.numV;

				double heDistance = getHeDistance(a);
				double heFactor = (double) (a - numHe) / dispersionHe;
				double vDistance = getVDistance(b);
				double vFactor = (double) (b - numV) / dispersionV;
				double pade = getPade(heDistance, vDistance);
				// A is the dissociating cluster
				dissPair.a00 += 1.0;
				dissPair.a01 += heFactor;
				dissPair.a02 += vFactor;
				dissPair.a10 += heDistance;
				dissPair.a11 += heDistance * heFactor;
				dissPair.a12 += heDistance * vFactor;
				dissPair.a20 += vDistance;
				dissPair.a21 += vDistance * heFactor;
				dissPair.a22 += vDistance * vFactor;
				dissPair.b0 += pade;
				dissPair.c0 += pade * heFactor;
				dissPair.d0 += pade * vFactor;
			});

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction,
		IReactant& disso) {

	// Check if we already know about the reaction.
	auto rkey = std::make_pair(&(reaction.first), &(reaction.second));
	auto it = effEmissionList.find(rkey);
	if (it == effEmissionList.end()) {

		// We did not already know about it.

		// Note that we emit from the two rectants according to the given
		// reaction.
		auto eret = effEmissionList.emplace(std::piecewise_construct,
				std::forward_as_tuple(rkey),
				std::forward_as_tuple(reaction,
						static_cast<PSICluster&>(reaction.first),
						static_cast<PSICluster&>(reaction.second)));
		// Since we already checked and didn't know about the reaction then,
		// we had better have added it with our emplace() call.
		assert(eret.second);
		it = eret.first;
	}
	assert(it != effEmissionList.end());
	auto& dissPair = it->second;

	// Values for grouping parameters
	int dissoLoHe = *(heBounds.begin()), dissoHiHe = *(heBounds.end()) - 1,
			dissoLoV = *(vBounds.begin()), dissoHiV = *(vBounds.end()) - 1,
			loHe = 0, hiHe = 0, loV = 0, hiV = 0, singleHeSize = 0,
			singleVSize = 0;

	if (dissPair.first.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(dissPair.first);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = dissPair.second.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}
	if (dissPair.second.getType() == ReactantType::PSISuper) {
		auto const& super = static_cast<PSISuperCluster const&>(dissPair.second);
		auto const& heBounds = super.getHeBounds();
		loHe = *(heBounds.begin());
		hiHe = *(heBounds.end()) - 1;
		auto const& vBounds = super.getVBounds();
		loV = *(vBounds.begin());
		hiV = *(vBounds.end()) - 1;
		auto singleComp = dissPair.first.getComposition();
		singleHeSize = singleComp[toCompIdx(Species::He)];
		singleVSize = singleComp[toCompIdx(Species::V)]
				- singleComp[toCompIdx(Species::I)]; // can be < 0
	}

	int heWidth = std::min(dissoHiHe, hiHe + singleHeSize)
			- std::max(dissoLoHe, loHe + singleHeSize) + 1;
	int vWidth = std::min(dissoHiV, hiV + singleVSize)
			- std::max(dissoLoV, loV + singleVSize) + 1;

	dissPair.a00 = heWidth * vWidth;

	dissPair.a01 = ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize), numHe);

	dissPair.a02 = ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize), numV);

	dissPair.a10 = ((double) (2 * vWidth) / (double) (sectionHeWidth - 1))
			* firstOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize), numHe);

	dissPair.a20 = ((double) (2 * heWidth) / (double) (sectionVWidth - 1))
			* firstOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize), numV);

	dissPair.a11 = ((double) (2 * vWidth)
			/ ((double) (sectionHeWidth - 1) * dispersionHe))
			* secondOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize), numHe);

	dissPair.a12 = ((double) (2 * vWidth) / (double) (sectionVWidth - 1))
			* firstOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize), numHe)
			* ((double) heWidth / dispersionV)
			* firstOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize), numV);

	dissPair.a21 = ((double) (2 * heWidth) / (double) (sectionHeWidth - 1))
			* firstOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize), numV)
			* ((double) vWidth / dispersionHe)
			* firstOrderSum(std::max(dissoLoHe, loHe + singleHeSize),
					std::min(dissoHiHe, hiHe + singleHeSize), numHe);

	dissPair.a22 = ((double) (2 * heWidth)
			/ ((double) (sectionVWidth - 1) * dispersionV))
			* secondOrderSum(std::max(dissoLoV, loV + singleVSize),
					std::min(dissoHiV, hiV + singleVSize), numV);

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
				/ ((double) nTot * (double) (sectionHeWidth - 1));

	if (sectionVWidth == 1)
		dispersionV = 1.0;
	else
		dispersionV = 2.0 * (nVSquare - (numV * (double) nTot * numV))
				/ ((double) nTot * (double) (sectionVWidth - 1));

	// Set the boundaries
	heBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>((numHe
					- (double) sectionHeWidth / 2.0) + 1),
			static_cast<IReactant::SizeType>((numHe
					- (double) sectionHeWidth / 2.0) + sectionHeWidth) + 1);
	vBounds = IntegerRange<IReactant::SizeType>(
			static_cast<IReactant::SizeType>((numV
					- (double) sectionVWidth / 2.0) + 1),
			static_cast<IReactant::SizeType>((numV
					- (double) sectionVWidth / 2.0) + sectionVWidth) + 1);

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

//			if (padeCoef.size() > 0) {
//				std::cout << name << " " << i << " " << j << " "
//						<< getPade(heDistance, vDistance) << std::endl;
//			}
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
				flux += value * (currPair.a00 * l0A + currPair.a10 * lHeA + currPair.a20 * lVA + currPair.b0);
				// Compute the momentum fluxes
				heMomentumFlux += value
				* (currPair.a01 * l0A + currPair.a11 * lHeA + currPair.a21 * lVA + currPair.c0);
				vMomentumFlux += value
				* (currPair.a02 * l0A + currPair.a12 * lHeA + currPair.a22 * lVA + currPair.d0);
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
				flux += value * (currPair.a00 * l0 + currPair.a10 * l1He + currPair.a20 * l1V + currPair.b0);
				// Compute the momentum fluxes
				heMomentumFlux -= value
				* (currPair.a01 * l0 + currPair.a11 * l1He + currPair.a21 * l1V + currPair.c0);
				vMomentumFlux -= value
				* (currPair.a02 * l0 + currPair.a12 * l1He + currPair.a22 * l1V + currPair.d0);
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
				* (currPair.a000 * l0A * l0B + currPair.a010 * l0A * lHeB
						+ currPair.a020 * l0A * lVB + currPair.a100 * lHeA * l0B
						+ currPair.a110 * lHeA * lHeB + currPair.a120 * lHeA * lVB
						+ currPair.a200 * lVA * l0B + currPair.a210 * lVA * lHeB
						+ currPair.a220 * lVA * lVB + currPair.b0 + currPair.b1 * l0A
						+ currPair.b2 * l0B + currPair.b3 * lHeA + currPair.b4 * lVA
						+ currPair.b5 * lHeB + currPair.b6 * lVB);
				// Compute the momentum fluxes
				heMomentumFlux += value
				* (currPair.a001 * l0A * l0B + currPair.a011 * l0A * lHeB
						+ currPair.a021 * l0A * lVB + currPair.a101 * lHeA * l0B
						+ currPair.a111 * lHeA * lHeB + currPair.a121 * lHeA * lVB
						+ currPair.a201 * lVA * l0B + currPair.a211 * lVA * lHeB
						+ currPair.a221 * lVA * lVB + currPair.c0 + currPair.c1 * l0A
						+ currPair.c2 * l0B + currPair.c3 * lHeA + currPair.c4 * lVA
						+ currPair.c5 * lHeB + currPair.c6 * lVB);
				vMomentumFlux += value
				* (currPair.a002 * l0A * l0B + currPair.a012 * l0A * lHeB
						+ currPair.a022 * l0A * lVB + currPair.a102 * lHeA * l0B
						+ currPair.a112 * lHeA * lHeB + currPair.a122 * lHeA * lVB
						+ currPair.a202 * lVA * l0B + currPair.a212 * lVA * lHeB
						+ currPair.a222 * lVA * lVB + currPair.d0 + currPair.d1 * l0A
						+ currPair.d2 * l0B + currPair.d3 * lHeA + currPair.d4 * lVA
						+ currPair.d5 * lHeB + currPair.d6 * lVB);
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
				* (currComb.a000 * l0B * l0 + currComb.a100 * l0B * l1He
						+ currComb.a200 * l0B * l1V + currComb.a010 * lHeB * l0
						+ currComb.a110 * lHeB * l1He + currComb.a210 * lHeB * l1V
						+ currComb.a020 * lVB * l0 + currComb.a120 * lVB * l1He
						+ currComb.a220 * lVB * l1V + currComb.b0 + currComb.b1 * l0
						+ currComb.b2 * l0B + currComb.b3 * l1He + currComb.b4 * l1V
						+ currComb.b5 * lHeB + currComb.b6 * lVB);
				// Compute the momentum fluxes
				heMomentumFlux -= value
				* (currComb.a001 * l0B * l0 + currComb.a101 * l0B * l1He
						+ currComb.a201 * l0B * l1V + currComb.a011 * lHeB * l0
						+ currComb.a111 * lHeB * l1He + currComb.a211 * lHeB * l1V
						+ currComb.a021 * lVB * l0 + currComb.a121 * lVB * l1He
						+ currComb.a221 * lVB * l1V + currComb.c0 + currComb.c1 * l0
						+ currComb.c2 * l0B + currComb.c3 * l1He + currComb.c4 * l1V
						+ currComb.c5 * lHeB + currComb.c6 * lVB);
				vMomentumFlux -= value
				* (currComb.a002 * l0B * l0 + currComb.a102 * l0B * l1He
						+ currComb.a202 * l0B * l1V + currComb.a012 * lHeB * l0
						+ currComb.a112 * lHeB * l1He + currComb.a212 * lHeB * l1V
						+ currComb.a022 * lVB * l0 + currComb.a122 * lVB * l1He
						+ currComb.a222 * lVB * l1V + currComb.d0 + currComb.d1 * l0
						+ currComb.d2 * l0B + currComb.d3 * l1He + currComb.d4 * l1V
						+ currComb.d5 * lHeB + currComb.d6 * lVB);
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
				* (currPair.a000 * l0B + currPair.a010 * lHeB + currPair.a020 * lVB + currPair.b1);
				heMomentumPartials[index] += value
				* (currPair.a001 * l0B + currPair.a011 * lHeB + currPair.a021 * lVB + currPair.c1);
				vMomentumPartials[index] += value
				* (currPair.a002 * l0B + currPair.a012 * lHeB + currPair.a022 * lVB + currPair.d1);
				index = firstReactant.getHeMomentumId() - 1;
				partials[index] += value
				* (currPair.a100 * l0B + currPair.a110 * lHeB + currPair.a120 * lVB + currPair.b3);
				heMomentumPartials[index] += value
				* (currPair.a101 * l0B + currPair.a111 * lHeB + currPair.a121 * lVB + currPair.c3);
				vMomentumPartials[index] += value
				* (currPair.a102 * l0B + currPair.a112 * lHeB + currPair.a122 * lVB + currPair.d3);
				index = firstReactant.getVMomentumId() - 1;
				partials[index] += value
				* (currPair.a200 * l0B + currPair.a210 * lHeB + currPair.a220 * lVB + currPair.b4);
				heMomentumPartials[index] += value
				* (currPair.a201 * l0B + currPair.a211 * lHeB + currPair.a221 * lVB + currPair.c4);
				vMomentumPartials[index] += value
				* (currPair.a202 * l0B + currPair.a212 * lHeB + currPair.a222 * lVB + currPair.d4);
				// Compute the contribution from the second part of the reacting pair
				index = secondReactant.getId() - 1;
				partials[index] += value
				* (currPair.a000 * l0A + currPair.a100 * lHeA + currPair.a200 * lVA + currPair.b2);
				heMomentumPartials[index] += value
				* (currPair.a001 * l0A + currPair.a101 * lHeA + currPair.a201 * lVA + currPair.c2);
				vMomentumPartials[index] += value
				* (currPair.a002 * l0A + currPair.a102 * lHeA + currPair.a202 * lVA + currPair.d2);
				index = secondReactant.getHeMomentumId() - 1;
				partials[index] += value
				* (currPair.a010 * l0A + currPair.a110 * lHeA + currPair.a210 * lVA + currPair.b5);
				heMomentumPartials[index] += value
				* (currPair.a011 * l0A + currPair.a111 * lHeA + currPair.a211 * lVA + currPair.c5);
				vMomentumPartials[index] += value
				* (currPair.a012 * l0A + currPair.a112 * lHeA + currPair.a212 * lVA + currPair.d5);
				index = secondReactant.getVMomentumId() - 1;
				partials[index] += value
				* (currPair.a020 * l0A + currPair.a120 * lHeA + currPair.a220 * lVA + currPair.b6);
				heMomentumPartials[index] += value
				* (currPair.a021 * l0A + currPair.a121 * lHeA + currPair.a221 * lVA + currPair.c6);
				vMomentumPartials[index] += value
				* (currPair.a022 * l0A + currPair.a122 * lHeA + currPair.a222 * lVA + currPair.d6);
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
				* (currComb.a000 * l0 + currComb.a100 * l1He + currComb.a200 * l1V + currComb.b2);
				heMomentumPartials[index] -= value
				* (currComb.a001 * l0 + currComb.a101 * l1He + currComb.a201 * l1V + currComb.c2);
				vMomentumPartials[index] -= value
				* (currComb.a002 * l0 + currComb.a102 * l1He + currComb.a202 * l1V + currComb.d2);
				index = cluster.getHeMomentumId() - 1;
				partials[index] -= value
				* (currComb.a010 * l0 + currComb.a110 * l1He + currComb.a210 * l1V + currComb.b5);
				heMomentumPartials[index] -= value
				* (currComb.a011 * l0 + currComb.a111 * l1He + currComb.a211 * l1V + currComb.c5);
				vMomentumPartials[index] -= value
				* (currComb.a012 * l0 + currComb.a112 * l1He + currComb.a212 * l1V + currComb.d5);
				index = cluster.getVMomentumId() - 1;
				partials[index] -= value
				* (currComb.a020 * l0 + currComb.a120 * l1He + currComb.a220 * l1V + currComb.b6);
				heMomentumPartials[index] -= value
				* (currComb.a021 * l0 + currComb.a121 * l1He + currComb.a221 * l1V + currComb.c6);
				vMomentumPartials[index] -= value
				* (currComb.a022 * l0 + currComb.a122 * l1He + currComb.a222 * l1V + currComb.d6);
				// Compute the contribution from this cluster
				index = id - 1;
				partials[index] -= value
				* (currComb.a000 * l0B + currComb.a010 * lHeB + currComb.a020 * lVB + currComb.b1);
				heMomentumPartials[index] -= value
				* (currComb.a001 * l0B + currComb.a011 * lHeB + currComb.a021 * lVB + currComb.c1);
				vMomentumPartials[index] -= value
				* (currComb.a002 * l0B + currComb.a012 * lHeB + currComb.a022 * lVB + currComb.d1);
				index = heMomId - 1;
				partials[index] -= value
				* (currComb.a100 * l0B + currComb.a110 * lHeB + currComb.a120 * lVB + currComb.b3);
				heMomentumPartials[index] -= value
				* (currComb.a101 * l0B + currComb.a111 * lHeB + currComb.a121 * lVB + currComb.c3);
				vMomentumPartials[index] -= value
				* (currComb.a102 * l0B + currComb.a112 * lHeB + currComb.a122 * lVB + currComb.d3);
				index = vMomId - 1;
				partials[index] -= value
				* (currComb.a200 * l0B + currComb.a210 * lHeB + currComb.a220 * lVB + currComb.b4);
				heMomentumPartials[index] -= value
				* (currComb.a201 * l0B + currComb.a211 * lHeB + currComb.a221 * lVB + currComb.c4);
				vMomentumPartials[index] -= value
				* (currComb.a202 * l0B + currComb.a212 * lHeB + currComb.a222 * lVB + currComb.d4);
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
				partials[index] += value * (currPair.a00);
				heMomentumPartials[index] += value * (currPair.a01);
				vMomentumPartials[index] += value * (currPair.a02);
				index = cluster.getHeMomentumId() - 1;
				partials[index] += value * (currPair.a10);
				heMomentumPartials[index] += value * (currPair.a11);
				vMomentumPartials[index] += value * (currPair.a12);
				index = cluster.getVMomentumId() - 1;
				partials[index] += value * (currPair.a20);
				heMomentumPartials[index] += value * (currPair.a21);
				vMomentumPartials[index] += value * (currPair.a22);
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
				partials[index] -= value * (currPair.a00);
				heMomentumPartials[index] -= value * (currPair.a01);
				vMomentumPartials[index] -= value * (currPair.a02);
				index = heMomId - 1;
				partials[index] -= value * (currPair.a10);
				heMomentumPartials[index] -= value * (currPair.a11);
				vMomentumPartials[index] -= value * (currPair.a12);
				index = vMomId - 1;
				partials[index] -= value * (currPair.a20);
				heMomentumPartials[index] -= value * (currPair.a21);
				vMomentumPartials[index] -= value * (currPair.a22);
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

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::ProductionCoefficientBase const& curr) const {

	os << "a[0-2][0-2][0-2]:" << ' ' << curr.a000 << ' ' << curr.a001 << ' '
			<< curr.a002 << ' ' << curr.a100 << ' ' << curr.a101 << ' '
			<< curr.a102 << ' ' << curr.a200 << ' ' << curr.a201 << ' '
			<< curr.a202 << ' ' << curr.a010 << ' ' << curr.a011 << ' '
			<< curr.a012 << ' ' << curr.a020 << ' ' << curr.a021 << ' '
			<< curr.a022 << ' ' << curr.a110 << ' ' << curr.a111 << ' '
			<< curr.a112 << ' ' << curr.a120 << ' ' << curr.a121 << ' '
			<< curr.a122 << ' ' << curr.a210 << ' ' << curr.a211 << ' '
			<< curr.a212 << ' ' << curr.a220 << ' ' << curr.a221 << ' '
			<< curr.a222;
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-2][0-2]:" << ' ' << currPair.a00 << ' ' << currPair.a01 << ' '
			<< currPair.a02 << ' ' << currPair.a10 << ' ' << currPair.a11 << ' '
			<< currPair.a12 << ' ' << currPair.a20 << ' ' << currPair.a21 << ' '
			<< currPair.a22;
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

