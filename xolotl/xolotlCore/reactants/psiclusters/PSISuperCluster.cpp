// Includes
#include <iterator>
#include "PSISuperCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <MathUtils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;


PSISuperCluster::PSISuperCluster(double _numHe, double _numV, int _nTot,
		int heWidth, int vWidth, IReactionNetwork& _network,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry, buildName(_numHe, _numV)), numHe(_numHe), numV(
				_numV), nTot(_nTot), heBounds(0, 0), vBounds(0, 0), l0(0.0), l1He(
				0.0), l1V(0.0), dispersionHe(0.0), dispersionV(0.0) {
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

	// Update the coefficients
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

	return;
}

void PSISuperCluster::participateIn(ProductionReaction& reaction, int a,
		int b) {

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
		const std::vector<PendingProductionReactionInfo>& pendingPRInfos) {

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

	// Update the coefficients
	std::for_each(pendingPRInfos.begin(), pendingPRInfos.end(),
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

	return;
}

void PSISuperCluster::participateIn(DissociationReaction& reaction, int a,
		int b, int c, int d) {

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

	// Update the coefficients
	double firstHeDistance = 0.0, firstVDistance = 0.0;
	if (reaction.dissociating.getType() == ReactantType::PSISuper) {
		auto const& super =
				static_cast<PSICluster const&>(reaction.dissociating);
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

	// Update the coefficients
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

	return;
}

void PSISuperCluster::emitFrom(DissociationReaction& reaction, int a, int b,
		int c, int d) {

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

	// Update the coeeficients
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

	return;
}

void PSISuperCluster::setHeVVector(const std::set<std::pair<int, int>>& vec) {
	// Copy the list of coordinates
	heVList = vec;

	// Initialize the dispersion sum
	double nHeSquare = 0.0, nVSquare = 0.0;
	// Update the network map, compute the radius and dispersions
    for (const auto& currPair : heVList) {
		double rad = (sqrt(3.0) / 4.0) * xolotlCore::tungstenLatticeConstant
				+ std::cbrt((3 * ipow<3>(xolotlCore::tungstenLatticeConstant) * currPair.second) / (8 * xolotlCore::pi))
				- std::cbrt((3 * ipow<3>(xolotlCore::tungstenLatticeConstant)) / (8 * xolotlCore::pi));
		reactionRadius += rad / (double) nTot;

		// Compute nSquare for the dispersion
		nHeSquare += (double) currPair.first * currPair.first;
		nVSquare += (double) currPair.second * currPair.second;
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
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, vDistance);
	}

	return conc;
}

double PSISuperCluster::getTotalHeliumConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster in the group times its number of helium
		conc += getConcentration(heDistance, vDistance) * (double) pair.first;
	}

	return conc;
}

double PSISuperCluster::getTotalVacancyConcentration() const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster in the group times its number of vacancies
		conc += getConcentration(heDistance, vDistance) * (double) pair.second;
	}

	return conc;
}

double PSISuperCluster::getIntegratedVConcentration(int v) const {
	// Initial declarations
	double heDistance = 0.0, vDistance = 0.0, conc = 0.0;

	// Loop on the indices
	for (auto const& pair : heVList) {
		// Skip the wrong V size
		if (pair.second != v) continue;

		// Compute the distances
		heDistance = getHeDistance(pair.first);
		vDistance = getVDistance(pair.second);

		// Add the concentration of each cluster
		conc += getConcentration(heDistance, vDistance);
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

	return;
}

PSISuperFlux PSISuperCluster::computeDissociationFlux() const {

	// Sum over all the dissociating pairs
    PSISuperFlux flux = 
        std::accumulate(effDissociatingList.begin(), effDissociatingList.end(),
            PSISuperFlux(),
			[this](const PSISuperFlux& running,
                        DissociationPairMap::value_type const& currMapItem) {

				auto const& currPair = currMapItem.second;

				// Get the dissociating clusters
				auto const& dissociatingCluster = currPair.first;
				double l0A = dissociatingCluster.getConcentration(0.0, 0.0);
				double lHeA = dissociatingCluster.getHeMomentum();
				double lVA = dissociatingCluster.getVMomentum();
				// Update the flux
				auto value = currPair.kConstant / (double) nTot;

                PSISuperFlux currFlux(
                    value * (currPair.a[0][0] * l0A + currPair.a[1][0] * lHeA + currPair.a[2][0] * lVA),
                    value * (currPair.a[0][1] * l0A + currPair.a[1][1] * lHeA + currPair.a[2][1] * lVA),
                    value * (currPair.a[0][2] * l0A + currPair.a[1][2] * lHeA + currPair.a[2][2] * lVA));

                return running + currFlux;
			});

	// Return the flux
	return flux;
}

PSISuperFlux PSISuperCluster::computeEmissionFlux() const {

	// Loop over all the emission pairs
	PSISuperFlux flux =
	    std::accumulate(effEmissionList.begin(), effEmissionList.end(),
            PSISuperFlux(),
			[this](const PSISuperFlux& running,
                        DissociationPairMap::value_type const& currMapItem) {

				auto const& currPair = currMapItem.second;

				// Update the flux
				auto value = currPair.kConstant / (double) nTot;

                PSISuperFlux currFlux(
                    value * (currPair.a[0][0] * l0 + currPair.a[1][0] * l1He + currPair.a[2][0] * l1V),
                    value * (currPair.a[0][1] * l0 + currPair.a[1][1] * l1He + currPair.a[2][1] * l1V),
                    value * (currPair.a[0][2] * l0 + currPair.a[1][2] * l1He + currPair.a[2][2] * l1V));

                return running + currFlux;
			});

	return flux;
}

PSISuperFlux PSISuperCluster::computeProductionFlux() const {

	// Sum over all the reacting pairs
	PSISuperFlux flux =
	    std::accumulate(effReactingList.begin(), effReactingList.end(),
            PSISuperFlux(),
			[this](const PSISuperFlux& running,
                        ProductionPairMap::value_type const& currMapItem) {

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

                PSISuperFlux currFlux(
                        value * (currPair.a[0][0][0] * l0A * l0B + currPair.a[0][1][0] * l0A * lHeB
						+ currPair.a[0][2][0] * l0A * lVB + currPair.a[1][0][0] * lHeA * l0B
						+ currPair.a[1][1][0] * lHeA * lHeB + currPair.a[1][2][0] * lHeA * lVB
						+ currPair.a[2][0][0] * lVA * l0B + currPair.a[2][1][0] * lVA * lHeB
						+ currPair.a[2][2][0] * lVA * lVB),
                value * (currPair.a[0][0][1] * l0A * l0B + currPair.a[0][1][1] * l0A * lHeB
						+ currPair.a[0][2][1] * l0A * lVB + currPair.a[1][0][1] * lHeA * l0B
						+ currPair.a[1][1][1] * lHeA * lHeB + currPair.a[1][2][1] * lHeA * lVB
						+ currPair.a[2][0][1] * lVA * l0B + currPair.a[2][1][1] * lVA * lHeB
						+ currPair.a[2][2][1] * lVA * lVB),
                value * (currPair.a[0][0][2] * l0A * l0B + currPair.a[0][1][2] * l0A * lHeB
						+ currPair.a[0][2][2] * l0A * lVB + currPair.a[1][0][2] * lHeA * l0B
						+ currPair.a[1][1][2] * lHeA * lHeB + currPair.a[1][2][2] * lHeA * lVB
						+ currPair.a[2][0][2] * lVA * l0B + currPair.a[2][1][2] * lVA * lHeB
						+ currPair.a[2][2][2] * lVA * lVB));

                return running + currFlux;
			});

	// Return the production flux
	return flux;
}

PSISuperFlux PSISuperCluster::computeCombinationFlux() const {

	// Sum over all the combining clusters
	PSISuperFlux flux =
	    std::accumulate(effCombiningList.begin(), effCombiningList.end(),
            PSISuperFlux(),
			[this](const PSISuperFlux& running,
                    CombiningClusterMap::value_type const& currMapItem) {

				// Get the combining cluster
				auto const& currComb = currMapItem.second;
				auto const& combiningCluster = currComb.first;
				double l0B = combiningCluster.getConcentration(0.0, 0.0);
				double lHeB = combiningCluster.getHeMomentum();
				double lVB = combiningCluster.getVMomentum();
				// Update the flux
				auto value = currComb.kConstant / (double) nTot;

                PSISuperFlux currFlux(
                    value * (currComb.a[0][0][0] * l0B * l0 + currComb.a[1][0][0] * l0B * l1He
						+ currComb.a[2][0][0] * l0B * l1V + currComb.a[0][1][0] * lHeB * l0
						+ currComb.a[1][1][0] * lHeB * l1He + currComb.a[2][1][0] * lHeB * l1V
						+ currComb.a[0][2][0] * lVB * l0 + currComb.a[1][2][0] * lVB * l1He
						+ currComb.a[2][2][0] * lVB * l1V),
                    value * (currComb.a[0][0][1] * l0B * l0 + currComb.a[1][0][1] * l0B * l1He
						+ currComb.a[2][0][1] * l0B * l1V + currComb.a[0][1][1] * lHeB * l0
						+ currComb.a[1][1][1] * lHeB * l1He + currComb.a[2][1][1] * lHeB * l1V
						+ currComb.a[0][2][1] * lVB * l0 + currComb.a[1][2][1] * lVB * l1He
						+ currComb.a[2][2][1] * lVB * l1V),
                    value * (currComb.a[0][0][2] * l0B * l0 + currComb.a[1][0][2] * l0B * l1He
						+ currComb.a[2][0][2] * l0B * l1V + currComb.a[0][1][2] * lHeB * l0
						+ currComb.a[1][1][2] * lHeB * l1He + currComb.a[2][1][2] * lHeB * l1V
						+ currComb.a[0][2][2] * lVB * l0 + currComb.a[1][2][2] * lVB * l1He
						+ currComb.a[2][2][2] * lVB * l1V));
                
                return running + currFlux;
			});

	return flux;
}


void PSISuperCluster::updateConcs(double* concs) const {

    // Compute our flux.
    auto flux = Reactant::computeFlux<PSISuperCluster>(*this);

    // Apply flux to current concentrations.
    updateConcsFromFlux(concs, flux);
}


void PSISuperCluster::updateConcsFromFlux(double* concs, const Flux& flux) const {
    // Have base class update concentrations from its part of flux.
    PSICluster::updateConcsFromFlux(concs, flux);

    // Update concentrations using our parts of the flux.
    auto const& superFlux = static_cast<PSISuperFlux const&>(flux);
    auto heIdx = getHeMomentumId() - 1;
    concs[heIdx] += superFlux.heMoment;
    auto vIdx = getVMomentumId() - 1;
    concs[vIdx] += superFlux.vMoment;
}


void PSISuperCluster::computePartialDerivatives(
		const std::vector<size_t>& startingIdx,
		const std::vector<int>& indices,
		std::vector<double>& vals) const {

    // Determine cluster's index into the size/indices/vals arrays.
    auto reactantIndex = getId() - 1;
    auto heReactantIndex = getHeMomentumId() - 1;
    auto vReactantIndex = getVMomentumId() - 1;

    // Get the inverse mappings from dense DOF space to
    // the indices/vals arrays.
	auto const& psiNetwork = static_cast<PSIClusterReactionNetwork const&>(network);
    auto const& partialsIdxMap = psiNetwork.getDFillInvMap(reactantIndex);
    auto const& hePartialsIdxMap = psiNetwork.getDFillInvMap(heReactantIndex);
    auto const& vPartialsIdxMap = psiNetwork.getDFillInvMap(vReactantIndex);

    // TODO do we want to wrap a vector around these?
    double* partials = &(vals[startingIdx[reactantIndex]]);
    double* hePartials = &(vals[startingIdx[heReactantIndex]]);
    double* vPartials = &(vals[startingIdx[vReactantIndex]]);

	// Get the partial derivatives for each reaction type
	computeProductionPartialDerivatives(partials, partialsIdxMap,
                                        hePartials, hePartialsIdxMap,
                                        vPartials, vPartialsIdxMap);
	computeCombinationPartialDerivatives(partials, partialsIdxMap,
                                        hePartials, hePartialsIdxMap,
                                        vPartials, vPartialsIdxMap);
	computeDissociationPartialDerivatives(partials, partialsIdxMap,
                                        hePartials, hePartialsIdxMap,
                                        vPartials, vPartialsIdxMap);
	computeEmissionPartialDerivatives(partials, partialsIdxMap,
                                        hePartials, hePartialsIdxMap,
                                        vPartials, vPartialsIdxMap);

	return;
}

void PSISuperCluster::computeProductionPartialDerivatives(
        double* partials,
        const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
        double* hePartials,
        const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
        double* vPartials,
        const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A

	// Loop over all the reacting pairs
	std::for_each(effReactingList.begin(), effReactingList.end(),
			[this,
                &partials, &partialsIdxMap,
                &hePartials, &hePartialsIdxMap,
                &vPartials, &vPartialsIdxMap](ProductionPairMap::value_type const& currMapItem) {

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
                {
                    auto index = firstReactant.getId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value
                    * (currPair.a[0][0][0] * l0B + currPair.a[0][1][0] * lHeB + currPair.a[0][2][0] * lVB);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value
                    * (currPair.a[0][0][1] * l0B + currPair.a[0][1][1] * lHeB + currPair.a[0][2][1] * lVB);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value
                    * (currPair.a[0][0][2] * l0B + currPair.a[0][1][2] * lHeB + currPair.a[0][2][2] * lVB);
                }
                {
                    auto index = firstReactant.getHeMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value
                    * (currPair.a[1][0][0] * l0B + currPair.a[1][1][0] * lHeB + currPair.a[1][2][0] * lVB);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value
                    * (currPair.a[1][0][1] * l0B + currPair.a[1][1][1] * lHeB + currPair.a[1][2][1] * lVB);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value
                    * (currPair.a[1][0][2] * l0B + currPair.a[1][1][2] * lHeB + currPair.a[1][2][2] * lVB);
                }
                {
                    auto index = firstReactant.getVMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value
                    * (currPair.a[2][0][0] * l0B + currPair.a[2][1][0] * lHeB + currPair.a[2][2][0] * lVB);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value
                    * (currPair.a[2][0][1] * l0B + currPair.a[2][1][1] * lHeB + currPair.a[2][2][1] * lVB);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value
                    * (currPair.a[2][0][2] * l0B + currPair.a[2][1][2] * lHeB + currPair.a[2][2][2] * lVB);
                }
                {
                    // Compute the contribution from the second part of the reacting pair
                    auto index = secondReactant.getId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value
                    * (currPair.a[0][0][0] * l0A + currPair.a[1][0][0] * lHeA + currPair.a[2][0][0] * lVA);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value
                    * (currPair.a[0][0][1] * l0A + currPair.a[1][0][1] * lHeA + currPair.a[2][0][1] * lVA);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value
                    * (currPair.a[0][0][2] * l0A + currPair.a[1][0][2] * lHeA + currPair.a[2][0][2] * lVA);
                }
                {
                    auto index = secondReactant.getHeMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value
                    * (currPair.a[0][1][0] * l0A + currPair.a[1][1][0] * lHeA + currPair.a[2][1][0] * lVA);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value
                    * (currPair.a[0][1][1] * l0A + currPair.a[1][1][1] * lHeA + currPair.a[2][1][1] * lVA);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value
                    * (currPair.a[0][1][2] * l0A + currPair.a[1][1][2] * lHeA + currPair.a[2][1][2] * lVA);
                }
                {
                    auto index = secondReactant.getVMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value
                    * (currPair.a[0][2][0] * l0A + currPair.a[1][2][0] * lHeA + currPair.a[2][2][0] * lVA);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value
                    * (currPair.a[0][2][1] * l0A + currPair.a[1][2][1] * lHeA + currPair.a[2][2][1] * lVA);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value
                    * (currPair.a[0][2][2] * l0A + currPair.a[1][2][2] * lHeA + currPair.a[2][2][2] * lVA);
                }
			});

	return;
}

void PSISuperCluster::computeCombinationPartialDerivatives(
        double* partials,
        const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
        double* hePartials,
        const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
        double* vPartials,
        const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A

	// Visit all the combining clusters
	std::for_each(effCombiningList.begin(), effCombiningList.end(),
			[this,
                &partials, &partialsIdxMap,
                &hePartials, &hePartialsIdxMap,
                &vPartials, &vPartialsIdxMap](CombiningClusterMap::value_type const& currMapItem) {
				// Get the combining clusters
				auto const& currComb = currMapItem.second;
				auto const& cluster = currComb.first;
				double l0B = cluster.getConcentration(0.0, 0.0);
				double lHeB = cluster.getHeMomentum();
				double lVB = cluster.getVMomentum();

				// Compute the contribution from the combining cluster
				auto value = currComb.kConstant / (double) nTot;
                {
                    auto index = cluster.getId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value
                    * (currComb.a[0][0][0] * l0 + currComb.a[1][0][0] * l1He + currComb.a[2][0][0] * l1V);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value
                    * (currComb.a[0][0][1] * l0 + currComb.a[1][0][1] * l1He + currComb.a[2][0][1] * l1V);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value
                    * (currComb.a[0][0][2] * l0 + currComb.a[1][0][2] * l1He + currComb.a[2][0][2] * l1V);
                }
                {
                    auto index = cluster.getHeMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value
                    * (currComb.a[0][1][0] * l0 + currComb.a[1][1][0] * l1He + currComb.a[2][1][0] * l1V);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value
                    * (currComb.a[0][1][1] * l0 + currComb.a[1][1][1] * l1He + currComb.a[2][1][1] * l1V);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value
                    * (currComb.a[0][1][2] * l0 + currComb.a[1][1][2] * l1He + currComb.a[2][1][2] * l1V);
                }
                {
                    auto index = cluster.getVMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value
                    * (currComb.a[0][2][0] * l0 + currComb.a[1][2][0] * l1He + currComb.a[2][2][0] * l1V);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value
                    * (currComb.a[0][2][1] * l0 + currComb.a[1][2][1] * l1He + currComb.a[2][2][1] * l1V);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value
                    * (currComb.a[0][2][2] * l0 + currComb.a[1][2][2] * l1He + currComb.a[2][2][2] * l1V);
                }
                {
                    // Compute the contribution from this cluster
                    auto index = id - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value
                    * (currComb.a[0][0][0] * l0B + currComb.a[0][1][0] * lHeB + currComb.a[0][2][0] * lVB);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value
                    * (currComb.a[0][0][1] * l0B + currComb.a[0][1][1] * lHeB + currComb.a[0][2][1] * lVB);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value
                    * (currComb.a[0][0][2] * l0B + currComb.a[0][1][2] * lHeB + currComb.a[0][2][2] * lVB);
                }
                {
                    auto index = heMomId - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value
                    * (currComb.a[1][0][0] * l0B + currComb.a[1][1][0] * lHeB + currComb.a[1][2][0] * lVB);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value
                    * (currComb.a[1][0][1] * l0B + currComb.a[1][1][1] * lHeB + currComb.a[1][2][1] * lVB);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value
                    * (currComb.a[1][0][2] * l0B + currComb.a[1][1][2] * lHeB + currComb.a[1][2][2] * lVB);
                }
                {
                    auto index = vMomId - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value
                    * (currComb.a[2][0][0] * l0B + currComb.a[2][1][0] * lHeB + currComb.a[2][2][0] * lVB);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value
                    * (currComb.a[2][0][1] * l0B + currComb.a[2][1][1] * lHeB + currComb.a[2][2][1] * lVB);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value
                    * (currComb.a[2][0][2] * l0B + currComb.a[2][1][2] * lHeB + currComb.a[2][2][2] * lVB);
                }
			});

	return;
}

void PSISuperCluster::computeDissociationPartialDerivatives(
        double* partials,
        const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
        double* hePartials,
        const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
        double* vPartials,
        const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)

	// Visit all the dissociating pairs
	std::for_each(effDissociatingList.begin(), effDissociatingList.end(),
			[this,
                &partials, &partialsIdxMap,
                &hePartials, &hePartialsIdxMap,
                &vPartials, &vPartialsIdxMap](DissociationPairMap::value_type const& currMapItem) {
				auto& currPair = currMapItem.second;

				// Get the dissociating clusters
				auto const& cluster = currPair.first;
				// Compute the contribution from the dissociating cluster
				auto value = currPair.kConstant / (double) nTot;
                {
                    auto index = cluster.getId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value * (currPair.a[0][0]);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value * (currPair.a[0][1]);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value * (currPair.a[0][2]);
                }
                {
                    auto index = cluster.getHeMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value * (currPair.a[1][0]);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value * (currPair.a[1][1]);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value * (currPair.a[1][2]);
                }
                {
                    auto index = cluster.getVMomentumId() - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] += value * (currPair.a[2][0]);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] += value * (currPair.a[2][1]);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] += value * (currPair.a[2][2]);
                }
			});

	return;
}

void PSISuperCluster::computeEmissionPartialDerivatives(
        double* partials,
        const ReactionNetwork::PartialsIdxMap& partialsIdxMap,
        double* hePartials,
        const ReactionNetwork::PartialsIdxMap& hePartialsIdxMap,
        double* vPartials,
        const ReactionNetwork::PartialsIdxMap& vPartialsIdxMap) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)

	// Visit all the emission pairs
	std::for_each(effEmissionList.begin(), effEmissionList.end(),
			[this,
                &partials, &partialsIdxMap,
                &hePartials, &hePartialsIdxMap,
                &vPartials, &vPartialsIdxMap](DissociationPairMap::value_type const& currMapItem) {
				auto& currPair = currMapItem.second;

				// Compute the contribution from the dissociating cluster
				auto value = currPair.kConstant / (double) nTot;
                {
                    auto index = id - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value * (currPair.a[0][0]);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value * (currPair.a[0][1]);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value * (currPair.a[0][2]);
                }
                {
                    auto index = heMomId - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value * (currPair.a[1][0]);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value * (currPair.a[1][1]);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value * (currPair.a[1][2]);
                }
                {
                    auto index = vMomId - 1;
                    auto partialsIdx = partialsIdxMap.at(index);
                    partials[partialsIdx] -= value * (currPair.a[2][0]);
                    auto hePartialsIdx = hePartialsIdxMap.at(index);
                    hePartials[hePartialsIdx] -= value * (currPair.a[2][1]);
                    auto vPartialsIdx = vPartialsIdxMap.at(index);
                    vPartials[vPartialsIdx] -= value * (currPair.a[2][2]);
                }
			});

	return;
}


void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::ProductionCoefficientBase const& curr) const {

	os << "a[0-2][0-2][0-2]: ";
    for(const auto& curr2D : curr.a) {
        for(const auto& curr1D : curr2D) {
            std::copy(curr1D.begin(), curr1D.end(),
                    std::ostream_iterator<double>(os, " "));
        }
    }
}

void PSISuperCluster::dumpCoefficients(std::ostream& os,
		PSISuperCluster::SuperClusterDissociationPair const& currPair) const {

	os << "a[0-2][0-2]: ";
    for(const auto& curr1D : currPair.a) {
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

