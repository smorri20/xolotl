#include <cassert>
#include <iterator>
#include "PSIClusterReactionNetwork.h"
#include "PSICluster.h"
#include "PSISuperCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

namespace xolotlCore {

PSIClusterReactionNetwork::PSIClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( {
                    ReactantType::V,
                    ReactantType::I,
                    ReactantType::He,
                    ReactantType::HeV,
                    ReactantType::HeI,
                    ReactantType::PSISuper
                },
                ReactantType::PSISuper, registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}


double PSIClusterReactionNetwork::calculateDissociationConstant(const DissociationReaction& reaction ) const {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// The atomic volume is computed by considering the BCC structure of the
	// tungsten. In a given lattice cell in tungsten there are tungsten atoms
	// at each corner and a tungsten atom in the center. The tungsten atoms at
	// the corners are shared across a total of eight cells. The fraction of
	// the volume of the lattice cell that is filled with tungsten atoms is the
	// atomic volume and is a_0^3/(8*1/8 + 1) = 0.5*a_0^3.
	double atomicVolume = 0.5 * xolotlCore::tungstenLatticeConstant
			* xolotlCore::tungstenLatticeConstant
			* xolotlCore::tungstenLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant;

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void PSIClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
    IReactant::SizeType firstSize = 0, secondSize = 0, productSize = 0;

	// Single species clustering (He, V, I)
	// We know here that only Xe_1 can cluster so we simplify the search
	// X_(a-i) + X_i --> X_a
	// Make a vector of types
    std::vector<ReactantType> typeVec { ReactantType::He, ReactantType::V, ReactantType::I };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {

        auto currType = *tvIter;

		// Get all the reactants of this type
		auto const& allTypeReactants = getAll(currType);
		// Loop on them
		for (auto firstIt = allTypeReactants.begin();
				firstIt != allTypeReactants.end(); firstIt++) {
			// Get its size
			firstSize = (*firstIt)->getSize();
			// Loop on the second cluster starting at the same pointer to avoid double counting
			for (auto secondIt = firstIt; secondIt != allTypeReactants.end();
					secondIt++) {
				// Get its size
				secondSize = (*secondIt)->getSize();
				productSize = firstSize + secondSize;
				// Get the product
				auto product = get(currType, productSize);
				// Check that the reaction can occur
				if (product
						&& ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt).get(), (*secondIt).get());
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
					product->createProduction(reaction);

					// Check if the reverse reaction is allowed
					checkDissociationConnectivity(product, reaction);
				}
			}
		}
	}

	// Helium absorption by HeV clusters
	// He_(a) + (He_b)(V_c) --> [He_(a+b)](V_c)
	// Get all the He and HeV clusters
	auto const& allHeReactants = getAll(ReactantType::He);
	auto const& allHeVReactants = getAll(ReactantType::HeV);
	auto const& allSuperReactants = getAll(ReactantType::PSISuper);
	// Loop on the He clusters
	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
			firstIt++) {
		// Skip if it can't diffuse
		if (xolotlCore::equal((*firstIt)->getDiffusionFactor(), 0.0))
			continue;
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allHeVReactants.begin();
				secondIt != allHeVReactants.end(); secondIt++) {
			// Get its composition
			auto& comp = (*secondIt)->getComposition();
			// Create the composition of the potential product
            auto newNumHe = comp[toCompIdx(Species::He)] + firstSize;
            auto newNumV = comp[toCompIdx(Species::V)];

			// Check if product already exists.
            IReactant::Composition newComp;
            newComp[toCompIdx(Species::He)] = newNumHe;
            newComp[toCompIdx(Species::V)] = newNumV;
			auto product = getCompound(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				// Check if it is a super cluster from the map
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
						(*secondIt).get());
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction, newNumHe, newNumV);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction,
                        newNumHe, newNumV);
			}
		}

		// Loop on the super clusters
		for (auto secondIt = allSuperReactants.begin();
				secondIt != allSuperReactants.end(); secondIt++) {
			auto superCluster = std::static_pointer_cast<PSISuperCluster>(*secondIt);
			IReactant * product = nullptr;
			// Get its boundaries
			auto boundaries = superCluster->getBoundaries();
			// Loop on them
			for (int i = boundaries[0]; i <= boundaries[1]; i++) {
				for (int j = boundaries[2]; j <= boundaries[3]; j++) {
					// Assume the product can only be a super cluster here
                    auto newNumHe = i + firstSize;
                    auto newNumV = j;
					product = getSuperFromComp(newNumHe, newNumV);
					// Check that the reaction can occur
					if (product
							&& ((*firstIt)->getDiffusionFactor() > 0.0
									|| (*secondIt)->getDiffusionFactor() > 0.0)) {
						// Create a production reaction
						auto reaction = std::make_shared<ProductionReaction>(
								(*firstIt).get(), (*secondIt).get());
						// Tell the reactants that they are in this reaction
						(*firstIt)->createCombination(reaction, i, j);
						(*secondIt)->createCombination(reaction, i, j);
						product->createProduction(reaction, newNumHe, newNumV,
                                i, j);

						// Check if the reverse reaction is allowed
						checkDissociationConnectivity(product, reaction,
								newNumHe, newNumV, i, j);
					}
				}
			}
		}
	}

	// Vacancy absorption by HeV clusters
	// (He_a)(V_b) + V_c --> (He_a)[V_(b+c)]
	// Get all the V clusters
	auto const& allVReactants = getAll(ReactantType::V);
	// Loop on the V clusters
	// Loop on the HeV clusters
	for (auto firstIt = allVReactants.begin(); firstIt != allVReactants.end();
			firstIt++) {
		// Skip if it can't diffuse
		if (xolotlCore::equal((*firstIt)->getDiffusionFactor(), 0.0))
			continue;
		// Get the V size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allHeVReactants.begin();
				secondIt != allHeVReactants.end(); secondIt++) {
			// Get its composition
			auto& comp = (*secondIt)->getComposition();
			// Create the composition of the potential product
            auto newNumHe = comp[toCompIdx(Species::He)];
            auto newNumV = comp[toCompIdx(Species::V)] + firstSize;

			// Check if product already exists.
            IReactant::Composition newComp;
            newComp[toCompIdx(Species::He)] = newNumHe;
            newComp[toCompIdx(Species::V)] = newNumV;
			auto product = getCompound(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
						(*secondIt).get());
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction, newNumHe, newNumV);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction,
						newNumHe, newNumV);
			}
		}

		// Loop on the super clusters
		for (auto secondIt = allSuperReactants.begin();
				secondIt != allSuperReactants.end(); secondIt++) {
			auto superCluster = std::static_pointer_cast<PSISuperCluster>(*secondIt);
			IReactant * product = nullptr;
			// Get its boundaries
			auto boundaries = superCluster->getBoundaries();
			// Loop on them
			for (int i = boundaries[0]; i <= boundaries[1]; i++) {
				for (int j = boundaries[2]; j <= boundaries[3]; j++) {

					// Assume the product can only be a super cluster here
                    auto newNumHe = i;
                    auto newNumV = j + firstSize;
					product = getSuperFromComp(newNumHe, newNumV);
					// Check that the reaction can occur
					if (product
							&& ((*firstIt)->getDiffusionFactor() > 0.0
									|| (*secondIt)->getDiffusionFactor() > 0.0)) {
						// Create a production reaction
						auto reaction = std::make_shared<ProductionReaction>(
								(*firstIt).get(), (*secondIt).get());
						// Tell the reactants that they are in this reaction
						(*firstIt)->createCombination(reaction, i, j);
						(*secondIt)->createCombination(reaction, i, j);
						product->createProduction(reaction,
                                newNumHe, newNumV,
                                i, j);

						// Check if the reverse reaction is allowed
						checkDissociationConnectivity(product, reaction,
                                newNumHe, newNumV,
                                i, j);
					}
				}
			}
		}
	}

	// Helium-Vacancy clustering
	// He_a + V_b --> (He_a)(V_b)
	// Loop on the He clusters
	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
			firstIt++) {
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allVReactants.begin();
				secondIt != allVReactants.end(); secondIt++) {
			// Get its size
			secondSize = (*secondIt)->getSize();
			// Create the composition of the potential product
            auto newNumHe = firstSize;
            auto newNumV = secondSize;

			// Get the product
            IReactant::Composition newComp;
            newComp[toCompIdx(Species::He)] = newNumHe;
            newComp[toCompIdx(Species::V)] = newNumV;
			auto product = getCompound(ReactantType::HeV, newComp);

			// Check if the product can be a super cluster
			if (!product) {
				product = getSuperFromComp(newNumHe, newNumV);
			}
			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
						(*secondIt).get());
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction, newNumHe, newNumV);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction,
						newNumHe, newNumV);
			}
		}
	}

	// Vacancy reduction by Interstitial absorption in HeV clusters
	// (He_a)(V_b) + (I_c) --> (He_a)[V_(b-c)]
	// Get all the I clusters
	auto const& allIReactants = getAll(ReactantType::I);
	// Loop on them
	for (auto firstIt = allIReactants.begin(); firstIt != allIReactants.end();
			firstIt++) {
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the HeV clusters
		for (auto secondIt = allHeVReactants.begin();
				secondIt != allHeVReactants.end(); secondIt++) {
			// Get its composition
			auto& comp = (*secondIt)->getComposition();
			// The product can be He or HeV
			IReactant * product = nullptr;
			if (comp[toCompIdx(Species::V)] == firstSize) {
				// The product is He
				product = get(ReactantType::He, comp[toCompIdx(Species::He)] );
			} else {
				// The product is HeV
				// Create the composition of the potential product
                IReactant::Composition newComp;
                newComp[toCompIdx(Species::He)] = comp[toCompIdx(Species::He)];
                newComp[toCompIdx(Species::V)] = comp[toCompIdx(Species::V)] - firstSize;
				// Get the product
				product = getCompound(ReactantType::HeV, newComp);
			}
			// Check that the reaction can occur
			if (product
					&& ((*firstIt)->getDiffusionFactor() > 0.0
							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
				// Create a production reaction
				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
						(*secondIt).get());
				// Tell the reactants that they are in this reaction
				(*firstIt)->createCombination(reaction);
				(*secondIt)->createCombination(reaction);
				product->createProduction(reaction);

				// Check if the reverse reaction is allowed
				checkDissociationConnectivity(product, reaction);
			}
		}

		// Loop on the super clusters
		for (auto secondIt = allSuperReactants.begin();
				secondIt != allSuperReactants.end(); secondIt++) {
			auto superCluster = std::static_pointer_cast<PSISuperCluster>(*secondIt);
			IReactant * product = nullptr;
			// Get its boundaries
			auto boundaries = superCluster->getBoundaries();
			// Loop on them
			for (auto i = boundaries[0]; i <= boundaries[1]; i++) {
				for (auto j = boundaries[2]; j <= boundaries[3]; j++) {
					// The product might be HeV or He
                    auto newNumHe = i;
                    auto newNumV = j - firstSize;

					// Get the product
					if (newNumV == 0) {
						// The product is He
						product = get(ReactantType::He, i);
					}
					else {
						// Create the composition of the potential product
                        IReactant::Composition newComp;
                        newComp[toCompIdx(Species::He)] = newNumHe;
                        newComp[toCompIdx(Species::V)] = newNumV;
						product = getCompound(ReactantType::HeV, newComp);

						// If the product doesn't exist check for super clusters
						if (!product) {
							product = getSuperFromComp(newNumHe, newNumV);
						}
					}
					// Check that the reaction can occur
					if (product
							&& ((*firstIt)->getDiffusionFactor() > 0.0
									|| (*secondIt)->getDiffusionFactor() > 0.0)) {
						// Create a production reaction
						auto reaction = std::make_shared<ProductionReaction>(
								(*firstIt).get(), (*secondIt).get());
						// Tell the reactants that they are in this reaction
						(*firstIt)->createCombination(reaction, i, j);
						(*secondIt)->createCombination(reaction, i, j);
						product->createProduction(reaction,
                                newNumHe, newNumV,
                                i, j);

						// Check if the reverse reaction is allowed
						checkDissociationConnectivity(product, reaction,
                                newNumHe, newNumV,
                                i, j);
					}
				}
			}
		}
	}

//	// Helium clustering leading to trap mutation
//	// He_a + He_b --> [He_(a+b)](V_c) + I_c
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the second He cluster starting at the same pointer to avoid double counting
//		for (auto secondIt = firstIt; secondIt != allHeReactants.end();
//				secondIt++) {
//			// Get its size
//			secondSize = (*secondIt)->getSize();
//			// Get the simple product
//			productSize = firstSize + secondSize;
//			auto product = get(ReactantType::He, productSize);
//			// Doesn't do anything if the product exist
//			if (product)
//				continue;
//
//			// Trap mutation is happening
//			// Loop on the possible I starting by the smallest
//			for (auto it = allIReactants.begin(); it != allIReactants.end();
//					it++) {
//				// Get the size of the I cluster
//				int iSize = (*it)->getSize();
//				// Create the composition of the potential product
//				std::vector<int> compositionVec = { firstSize + secondSize,
//						iSize, 0 };
//				product = getCompound(ReactantType::HeV, compositionVec);
//				// Check that the reaction can occur
//				if (product
//						&& ((*firstIt)->getDiffusionFactor() > 0.0
//								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//					// Create a production reaction
//					auto reaction = std::make_shared<ProductionReaction>(
//							(*firstIt).get(), (*secondIt).get());
//					// Tell the reactants that they are in this reaction
//					(*firstIt)->createCombination(reaction);
//					(*secondIt)->createCombination(reaction);
//					product->createProduction(reaction);
//					(*it)->createProduction(reaction);
//
//					// Stop the loop on I clusters here
//					break;
//				}
//			}
//		}
//	}

//	// Helium absorption by HeV leading to trap mutation
//	// (He_a)(V_b) + He_c --> [He_(a+c)][V_(b+d)] + I_d
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the HeV clusters
//		for (auto secondIt = allHeVReactants.begin();
//				secondIt != allHeVReactants.end(); secondIt++) {
//			// Get its composition
//			auto& comp = (*secondIt)->getComposition();
//			// Get the simple product
//			std::vector<int> compositionVec = { firstSize + comp[toCompIdx(Species::He)],
//					comp[toCompIdx(Species::V)], 0 };
//			auto product = getCompound(ReactantType::HeV, compositionVec);
//			// Doesn't do anything if the product exist
//			if (product)
//				continue;
//
//			// Trap mutation is happening
//			// Loop on the possible I starting by the smallest
//			for (auto it = allIReactants.begin(); it != allIReactants.end();
//					it++) {
//				// Get the size of the I cluster
//				int iSize = (*it)->getSize();
//				// Create the composition of the potential product
//				compositionVec[1] = comp[toCompIdx(Species::V)] + iSize;
//				product = getCompound(ReactantType::HeV, compositionVec);
//				// Check that the reaction can occur
//				if (product
//						&& ((*firstIt)->getDiffusionFactor() > 0.0
//								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//					// Create a production reaction
//					auto reaction = std::make_shared<ProductionReaction>(
//							(*firstIt).get(), (*secondIt).get());
//					// Tell the reactants that they are in this reaction
//					(*firstIt)->createCombination(reaction);
//					(*secondIt)->createCombination(reaction);
//					product->createProduction(reaction);
//					(*it)->createProduction(reaction);
//
//					// Stop the loop on I clusters here
//					break;
//				}
//			}
//		}
//	}

	// Vacancy-Interstitial annihilation
	// I_a + V_b
	//        --> I_(a-b), if a > b
	//        --> V_(b-a), if a < b
	//        --> 0, if a = b
	// Loop on the I clusters
	for (auto firstIt = allIReactants.begin(); firstIt != allIReactants.end();
			firstIt++) {
		// Get its size
		firstSize = (*firstIt)->getSize();
		// Loop on the V clusters
		for (auto secondIt = allVReactants.begin();
				secondIt != allVReactants.end(); secondIt++) {
			// Get its size
			secondSize = (*secondIt)->getSize();
			// Check the possibilities
			if (firstSize > secondSize) {
				// Get the product
				productSize = firstSize - secondSize;
				auto product = get(ReactantType::I, productSize);
				// Check that the reaction can occur
				if (product
						&& ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt).get(), (*secondIt).get());
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
					product->createProduction(reaction);
				}
			} else if (firstSize < secondSize) {
				// Get the product
				productSize = secondSize - firstSize;
				auto product = get(ReactantType::V, productSize);
				// Check that the reaction can occur
				if (product
						&& ((*firstIt)->getDiffusionFactor() > 0.0
								|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt).get(), (*secondIt).get());
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
					product->createProduction(reaction);
				}

			} else {
				// Annihilation
				// Check that the reaction can occur
				if (((*firstIt)->getDiffusionFactor() > 0.0
						|| (*secondIt)->getDiffusionFactor() > 0.0)) {
					// Create a production reaction
					auto reaction = std::make_shared<ProductionReaction>(
							(*firstIt).get(), (*secondIt).get());
					// Tell the reactants that they are in this reaction
					(*firstIt)->createCombination(reaction);
					(*secondIt)->createCombination(reaction);
				}
			}
		}
	}

//	// Helium absorption by HeI clusters
//	// He_(a) + (He_b)(I_c) --> [He_(a+b)](I_c)
//	// Get all the HeI clusters
//	auto const& allHeIReactants = getAll(ReactantType::HeI);
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the HeV clusters
//		for (auto secondIt = allHeIReactants.begin();
//				secondIt != allHeIReactants.end(); secondIt++) {
//			// Get its composition
//			auto& comp = (*secondIt)->getComposition();
//			// Create the composition of the potential product
//			std::vector<int> compositionVec = { comp[toCompIdx(Species::He)] + firstSize, 0,
//					comp[toCompIdx(Species::I)] };
//			// Get the product
//			auto product = getCompound(ReactantType::HeI, compositionVec);
//			// Check that the reaction can occur
//			if (product
//					&& ((*firstIt)->getDiffusionFactor() > 0.0
//							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//				// Create a production reaction
//				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
//						(*secondIt).get());
//				// Tell the reactants that they are in this reaction
//				(*firstIt)->createCombination(reaction);
//				(*secondIt)->createCombination(reaction);
//				product->createProduction(reaction);
//
//				// Check if the reverse reaction is allowed
//				checkDissociationConnectivity(product, reaction);
//			}
//		}
//	}
//
//	// Single Interstitial absorption by HeI clusters
//	// (He_a)(I_b) + I --> (He_a)[I_(b+1)]
//	// Get the single interstitial cluster
//	auto singleInterstitialCluster = get(ReactantType::I, 1);
//	// Loop on the HeI clusters
//	for (auto secondIt = allHeIReactants.begin();
//			secondIt != allHeIReactants.end(); secondIt++) {
//		// Get its composition
//		auto& comp = (*secondIt)->getComposition();
//		// Create the composition of the potential product
//		std::vector<int> compositionVec = { comp[toCompIdx(Species::He)], 0, comp[toCompIdx(Species::I)] + 1 };
//		// Get the product
//		auto product = getCompound(ReactantType::HeI, compositionVec);
//		// Check that the reaction can occur
//		if (product
//				&& (singleInterstitialCluster->getDiffusionFactor() > 0.0
//						|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//			// Create a production reaction
//			auto reaction = std::make_shared<ProductionReaction>(
//					singleInterstitialCluster, (*secondIt).get());
//			// Tell the reactants that they are in this reaction
//			singleInterstitialCluster->createCombination(reaction);
//			(*secondIt)->createCombination(reaction);
//			product->createProduction(reaction);
//
//			// Check if the reverse reaction is allowed
//			checkDissociationConnectivity(product, reaction);
//		}
//	}
//
//	// Helium-Interstitial clustering
//	// He_a + I_b --> (He_a)(I_b)
//	// Loop on the He clusters
//	for (auto firstIt = allHeReactants.begin(); firstIt != allHeReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the I clusters
//		for (auto secondIt = allIReactants.begin();
//				secondIt != allIReactants.end(); secondIt++) {
//			// Get its size
//			secondSize = (*secondIt)->getSize();
//			// Create the composition of the potential product
//			std::vector<int> compositionVec = { firstSize, 0, secondSize };
//			// Get the product
//			auto product = getCompound(ReactantType::HeI, compositionVec);
//			// Check that the reaction can occur
//			if (product
//					&& ((*firstIt)->getDiffusionFactor() > 0.0
//							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//				// Create a production reaction
//				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
//						(*secondIt).get());
//				// Tell the reactants that they are in this reaction
//				(*firstIt)->createCombination(reaction);
//				(*secondIt)->createCombination(reaction);
//				product->createProduction(reaction);
//
//				// Check if the reverse reaction is allowed
//				checkDissociationConnectivity(product, reaction);
//			}
//		}
//	}
//
//	// Interstitial reduction by Vacancy absorption in HeI clusters
//	// (He_a)(I_b) + (V_c) --> (He_a)[I_(b-c)]
//	// Loop on V clusters
//	for (auto firstIt = allVReactants.begin(); firstIt != allVReactants.end();
//			firstIt++) {
//		// Get its size
//		firstSize = (*firstIt)->getSize();
//		// Loop on the HeI clusters
//		for (auto secondIt = allHeIReactants.begin();
//				secondIt != allHeIReactants.end(); secondIt++) {
//			// Get its composition
//			auto& comp = (*secondIt)->getComposition();
//			// The product can be He or HeI
//			IReactant * product = nullptr;
//			if (comp[toCompIdx(Species::I)] == firstSize) {
//				// The product is He
//				product = get(ReactantType::He, comp[toCompIdx(Species::He)] );
//			} else {
//				// The product is HeI
//				// Create the composition of the potential product
//				std::vector<int> compositionVec = { comp[toCompIdx(Species::He)] , 0, comp[toCompIdx(Species::I)]
//						- firstSize };
//				// Get the product
//				product = getCompound(ReactantType::HeI, compositionVec);
//			}
//			// Check that the reaction can occur
//			if (product
//					&& ((*firstIt)->getDiffusionFactor() > 0.0
//							|| (*secondIt)->getDiffusionFactor() > 0.0)) {
//				// Create a production reaction
//				auto reaction = std::make_shared<ProductionReaction>((*firstIt).get(),
//						(*secondIt));
//				// Tell the reactants that they are in this reaction
//				(*firstIt)->createCombination(reaction);
//				(*secondIt)->createCombination(reaction);
//				product->createProduction(reaction);
//
//				// Check if the reverse reaction is allowed
//				checkDissociationConnectivity(product, reaction);
//			}
//		}
//	}

	return;
}

void PSIClusterReactionNetwork::checkDissociationConnectivity(
		IReactant * emittingReactant,
		std::shared_ptr<ProductionReaction> reaction, int a, int b, int c,
		int d) {
	// Check if at least one of the potentially emitted cluster is size one
	if (reaction->first->getSize() != 1 && reaction->second->getSize() != 1) {
		// Don't add the reverse reaction
		return;
	}
	// remove He+He
	if (reaction->first->getSize() == 1 && reaction->second->getSize() == 1
			&& reaction->first->getType() == ReactantType::He
			&& reaction->second->getType() == ReactantType::He) {
		// Don't add the reverse reaction
		return;
	}

//	// Check for trap mutations (with XOR)
//	if ((reaction->first->getType() == ReactantType::I)
//			== !(reaction->second->getType() == ReactantType::I)) {
//		// Don't add the reverse reaction
//		return;
//	}

	// The reaction can occur, create the dissociation
	// Create a dissociation reaction
	auto dissociationReaction = std::make_shared<DissociationReaction>(
			emittingReactant, reaction->first, reaction->second);
	// Set the reverse reaction
	dissociationReaction->reverseReaction = reaction.get();
	// Tell the reactants that their are in this reaction
	reaction->first->createDissociation(dissociationReaction, a, b, c, d);
	reaction->second->createDissociation(dissociationReaction, a, b, c, d);
	emittingReactant->createEmission(dissociationReaction, a, b, c, d);

	return;
}

void PSIClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	computeRateConstants();

	return;
}

double PSIClusterReactionNetwork::getTemperature() const {
	return temperature;
}

IReactant * PSIClusterReactionNetwork::get(ReactantType type,
		IReactant::SizeType size) const {

	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	// Only pull the reactant if the name and size are valid
	if ((type == ReactantType::He || type == ReactantType::V || type == ReactantType::I) && size >= 1) {

        IReactant::Composition composition;
        composition[toCompIdx(toSpecies(type))] = size;

		//std::string encodedName = PSICluster::encodeCompositionAsName(composition);
		// Check if the reactant is in the map
        auto iter = singleSpeciesMap.find(composition);
        if (iter != singleSpeciesMap.end()) {
            retReactant = iter->second;
		}
	}

	return retReactant.get();
}

IReactant * PSIClusterReactionNetwork::getCompound(ReactantType type,
        const IReactant::Composition& comp) const {

	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if (type == ReactantType::HeV || type == ReactantType::HeI) {

		// Check if the reactant is in the map
        auto iter = mixedSpeciesMap.find(comp);
        if (iter != mixedSpeciesMap.end()) {
			retReactant = iter->second;
		}
	}

	return retReactant.get();
}

IReactant * PSIClusterReactionNetwork::getSuper(ReactantType type,
        const IReactant::Composition& comp) const {

	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if (type == ReactantType::PSISuper) {

		// Check if the reactant is in the map
        auto iter = superSpeciesMap.find(comp);
        if (iter != superSpeciesMap.end()) {
			retReactant = iter->second;
		}
	}

	return retReactant.get();
}


void PSIClusterReactionNetwork::add(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto& composition = reactant->getComposition();

		// Get the species sizes
		numHe = composition[toCompIdx(Species::He)] ;
		numV = composition[toCompIdx(Species::V)] ;
		numI = composition[toCompIdx(Species::I)] ;

		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
		if (isMixed && mixedSpeciesMap.count(composition) == 0) {
			// Put the compound in its map
            mixedSpeciesMap.emplace(composition, reactant);
		} else if (!isMixed && singleSpeciesMap.count(composition) == 0) {
			/// Put the reactant in its map
            singleSpeciesMap.emplace(composition, reactant);
		} else {
			std::stringstream errStream;
			errStream << "PSIClusterReactionNetwork Message: "
					<< "Duplicate Reactant (He=" << numHe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the max cluster size for the new cluster's type.
        IReactant::SizeType clusterSize = numHe + numV + numI;
        if(clusterSize > maxClusterSizeMap[reactant->getType()]) {
            maxClusterSizeMap[reactant->getType()] = clusterSize;
        }

		// Set the id for this cluster
        // (It is networkSize+1 because we haven't added it to the network yet.)
		reactant->setId(size()+1);

        // Add reactant to our per-type map.
        clusterTypeMap.at(reactant->getType()).push_back(reactant);

		// Add the pointer to the list of all clusters
        allReactants.push_back(reactant.get());
	}

	return;
}

void PSIClusterReactionNetwork::addSuper(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numHe = 0, numV = 0, numI = 0;
	bool isMixed = false;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto& composition = reactant->getComposition();

		// Get the species sizes
		numHe = composition[toCompIdx(Species::He)] ;
		numV = composition[toCompIdx(Species::V)] ;
		numI = composition[toCompIdx(Species::I)] ;
		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numHe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
		if (isMixed && superSpeciesMap.count(composition) == 0) {
			// Put the compound in its map
            superSpeciesMap.emplace(composition, reactant);
		} else {
			std::stringstream errStream;
			errStream << "PSIClusterReactionNetwork Message: "
					<< "Duplicate Super Reactant (He=" << numHe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}

		// Set the id for this cluster
        // (It is networkSize+1 because we haven't added it to the network yet.)
		reactant->setId(size() + 1);
		// Add cluster to our per-type map.
        clusterTypeMap.at(reactant->getType()).push_back(reactant);

		// Add the pointer to the list of all clusters
		allReactants.push_back(reactant.get());
	}

	return;
}


void PSIClusterReactionNetwork::buildSuperClusterIndex(const std::vector<IReactant::SizeType>& bounds) {

    // Save the bounds to use.
    boundVector = bounds;

    // Build a map of super clusters, keyed by (baseHe, baseV) pairs
    // where base* indicates the lower bound of the super cluster's
    // interval for that species type.
    auto const& superClusters = clusterTypeMap[ReactantType::PSISuper];
    std::for_each(superClusters.begin(), superClusters.end(),
        [this](const std::shared_ptr<IReactant>& currCluster) {
            // Add the super cluster to our lookup map based on
            // its He and V intervals.
            auto const& comp = currCluster->getComposition();
            auto currHe = comp[toCompIdx(Species::He)];
            auto currV = comp[toCompIdx(Species::V)];
            auto heBoundsIntervalBase = findBoundsIntervalBase(currHe);
            auto vBoundsIntervalBase = findBoundsIntervalBase(currV);

            superClusterLookupMap.emplace(std::make_pair(heBoundsIntervalBase, vBoundsIntervalBase), *currCluster);
        });
}


void PSIClusterReactionNetwork::removeReactants(
		const IReactionNetwork::ReactantVector& doomedReactants) {

	// Build a ReactantMatcher functor for the doomed reactants.
	// Doing this here allows us to construct the canonical composition
	// strings for the doomed reactants once and reuse them.
	// If we used an anonymous functor object in the std::remove_if
	// calls we would build these strings several times in this function.
	ReactionNetwork::ReactantMatcher doomedReactantMatcher(doomedReactants);

	// Remove the doomed reactants from our collection of all known reactants.
	auto ariter = std::remove_if(allReactants.begin(), allReactants.end(),
			doomedReactantMatcher);
	allReactants.erase(ariter, allReactants.end());

	// Remove the doomed reactants from the type-specific cluster vectors.
	// First, determine all cluster types used by clusters in the collection
	// of doomed reactants...
	std::set<ReactantType> typesUsed;
	for (auto reactant : doomedReactants) {
		typesUsed.insert(reactant->getType());
	}

	// ...Next, examine each type's collection of clusters and remove the
	// doomed reactants.
	for (auto currType : typesUsed) {
		auto& clusters = clusterTypeMap[currType];
		auto citer = std::remove_if(clusters.begin(), clusters.end(),
				doomedReactantMatcher);
		clusters.erase(citer, clusters.end());
	}

	// Remove the doomed reactants from the SpeciesMap.
	// We cannot use std::remove_if and our ReactantMatcher here
	// because std::remove_if reorders the elements in the underlying
	// container to move the doomed elements to the end of the container,
	// but the std::map doesn't support reordering.
	for (auto reactant : doomedReactants) {
		if (reactant->isMixed()) {
			mixedSpeciesMap.erase(reactant->getComposition());
        }
		else {
			singleSpeciesMap.erase(reactant->getComposition());
        }
	}

	return;
}

void PSIClusterReactionNetwork::reinitializeNetwork() {

	// Reset the Ids
	int id = 0;
	for (auto it = allReactants.begin(); it != allReactants.end(); ++it) {
		id++;
		(*it)->setId(id);
		(*it)->setHeMomentumId(id);
		(*it)->setVMomentumId(id);
	}

	// Get all the super clusters and loop on them
    for (auto& currCluster : clusterTypeMap[ReactantType::PSISuper]) {

		id++;
		currCluster->setHeMomentumId(id);
		id++;
		currCluster->setVMomentumId(id);

		// Update the HeV size
		auto cluster = (PSISuperCluster *) currCluster.get();
		auto bounds = cluster->getBoundaries();
        IReactant::SizeType clusterSize = bounds[1] + bounds[3];
        if(clusterSize > maxClusterSizeMap[ReactantType::HeV]) {
            maxClusterSizeMap[ReactantType::HeV] = clusterSize;
        }
	}

	return;
}

void PSIClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	for (auto it = allReactants.begin(); it != allReactants.end(); ++it) {
		(*it)->resetConnectivities();
	}

	return;
}

void PSIClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {
	// Local Declarations
	int size = allReactants.size();
	int id = 0;

	// Set the concentrations
	concUpdateCounter->increment();	// increment the update concentration counter
	for (int i = 0; i < size; i++) {
		id = allReactants.at(i)->getId() - 1;
		allReactants.at(i)->setConcentration(concentrations[id]);
	}

	// Set the moments
	for (int i = size - getAll(ReactantType::PSISuper).size(); i < size; i++) {
		auto cluster = (PSISuperCluster *) allReactants.at(i);
		id = cluster->getId() - 1;
		cluster->setZerothMomentum(concentrations[id]);
		id = cluster->getHeMomentumId() - 1;
		cluster->setHeMomentum(concentrations[id]);
		id = cluster->getVMomentumId() - 1;
		cluster->setVMomentum(concentrations[id]);
	}

	return;
}

void PSIClusterReactionNetwork::getDiagonalFill(int *diagFill) {
	// Get all the super clusters
	auto const& superClusters = getAll(ReactantType::PSISuper);

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Declarations for the loop
	std::vector<int> connectivity;
	int connectivityLength, id, index;

	// Get the connectivity for each reactant
	for (int i = 0; i < size(); i++) {
		// Get the reactant and its connectivity
		auto& reactant = allReactants.at(i);
		connectivity = reactant->getConnectivity();
		connectivityLength = connectivity.size();
		// Get the reactant id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getId() - 1;
		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = id * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}
	// Get the connectivity for each moment
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the reactant and its connectivity
		auto reactant = superClusters[i];
		connectivity = reactant->getConnectivity();
		connectivityLength = connectivity.size();
		// Get the helium momentum id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getHeMomentumId() - 1;

		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = (id) * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;

		// Get the vacancy momentum id so that the connectivity can be lined up in
		// the proper column
		id = reactant->getVMomentumId() - 1;

		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			index = (id) * dof + j;
			diagFill[index] = connectivity[j];
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

double PSIClusterReactionNetwork::getTotalAtomConcentration() {
	// Initial declarations
	double heliumConc = 0.0;

	// Get all the He clusters
	auto const& heClusters = getAll(ReactantType::He);
	// Loop on them
	for (int i = 0; i < heClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heClusters[i];
		double size = cluster->getSize();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster->getConcentration() * size;
	}

	// Get all the HeV clusters
	auto const& heVClusters = getAll(ReactantType::HeV);
	// Loop on them
	for (int i = 0; i < heVClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heVClusters[i];
		auto& comp = cluster->getComposition();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster->getConcentration() * comp[toCompIdx(Species::He)] ;
	}

	// Get all the super clusters
	auto const& superClusters = getAll(ReactantType::PSISuper);
	// Loop on them
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the cluster
		auto cluster = std::static_pointer_cast<PSISuperCluster>(superClusters[i]);

		// Add its total helium concentration helium concentration
		heliumConc += cluster->getTotalHeliumConcentration();
	}

	return heliumConc;
}

double PSIClusterReactionNetwork::getTotalTrappedAtomConcentration() {
	// Initial declarations
	double heliumConc = 0.0;

	// Get all the HeV clusters
	auto const& heVClusters = getAll(ReactantType::HeV);
	// Loop on them
	for (int i = 0; i < heVClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heVClusters[i];
		auto& comp = cluster->getComposition();

		// Add the concentration times the He content to the total helium concentration
		heliumConc += cluster->getConcentration() * comp[toCompIdx(Species::He)] ;
	}

	// Get all the super clusters
	auto const& superClusters = getAll(ReactantType::PSISuper);
	// Loop on them
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the cluster
		auto cluster = std::static_pointer_cast<PSISuperCluster>(superClusters[i]);

		// Add its total helium concentration
		heliumConc += cluster->getTotalHeliumConcentration();
	}

	return heliumConc;
}

double PSIClusterReactionNetwork::getTotalVConcentration() {
	// Initial declarations
	double vConc = 0.0;

	// Get all the V clusters
	auto const& vClusters = getAll(ReactantType::V);
	// Loop on them
	for (int i = 0; i < vClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = vClusters[i];
		double size = cluster->getSize();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster->getConcentration() * size;
	}

	// Get all the HeV clusters
	auto const& heVClusters = getAll(ReactantType::HeV);
	// Loop on them
	for (int i = 0; i < heVClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = heVClusters[i];
		auto& comp = cluster->getComposition();

		// Add the concentration times the V content to the total vacancy concentration
		vConc += cluster->getConcentration() * comp[toCompIdx(Species::V)] ;
	}

	// Get all the super clusters
	auto const& superClusters = getAll(ReactantType::PSISuper);
	// Loop on them
	for (int i = 0; i < superClusters.size(); i++) {
		// Get the cluster
		auto cluster = std::static_pointer_cast<PSISuperCluster>(superClusters[i]);

		// Add its total vacancy concentration
		vConc += cluster->getTotalVacancyConcentration();
	}

	return vConc;
}

double PSIClusterReactionNetwork::getTotalIConcentration() {
	// Initial declarations
	double iConc = 0.0;

	// Get all the V clusters
	auto const& iClusters = getAll(ReactantType::I);
	// Loop on them
	for (int i = 0; i < iClusters.size(); i++) {
		// Get the cluster and its composition
		auto cluster = iClusters[i];
		double size = cluster->getSize();

		// Add the concentration times the I content to the total interstitial concentration
		iConc += cluster->getConcentration() * size;
	}

	return iConc;
}

void PSIClusterReactionNetwork::computeRateConstants() {
	// Local declarations
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on all the production reactions
    for (auto& currReactionInfo : productionReactionMap) {

        auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateReactionRateConstant(*currReaction);
		// Set it in the reaction
		currReaction->kConstant = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on all the dissociation reactions
    for (auto& currReactionInfo : dissociationReactionMap) {

        auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateDissociationConstant(*currReaction);

		// Set it in the reaction
		currReaction->kConstant = rate;
	}

	// Set the biggest rate
	biggestRate = biggestProductionRate;

	return;
}

void PSIClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset) {
	// Initial declarations
	double flux = 0.0;
	int reactantIndex = 0;
	auto const& superClusters = getAll(ReactantType::PSISuper);

	// ----- Compute all of the new fluxes -----
	for (int i = 0; i < size(); i++) {
		auto& cluster = allReactants.at(i);
		// Compute the flux
		flux = cluster->getTotalFlux();
		// Update the concentration of the cluster
		reactantIndex = cluster->getId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	// ---- Moments ----
	for (int i = 0; i < superClusters.size(); i++) {
		auto superCluster = std::static_pointer_cast<PSISuperCluster>(superClusters[i]);

		// Compute the helium momentum flux
		flux = superCluster->getHeMomentumFlux();
		// Update the concentration of the cluster
		reactantIndex = superCluster->getHeMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;

		// Compute the vacancy momentum flux
		flux = superCluster->getVMomentumFlux();
		// Update the concentration of the cluster
		reactantIndex = superCluster->getVMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	return;
}

void PSIClusterReactionNetwork::computeAllPartials(double *vals, int *indices,
		int *size) {
	// Initial declarations
	int reactantIndex = 0, pdColIdsVectorSize = 0;
	const int dof = getDOF();
	std::vector<double> clusterPartials;
	clusterPartials.resize(dof, 0.0);
	// Get the super clusters
	auto const& superClusters = getAll(ReactantType::PSISuper);

	// Update the column in the Jacobian that represents each normal reactant
	for (int i = 0; i < this->size() - superClusters.size(); i++) {
		auto& reactant = allReactants.at(i);
		// Get the reactant index
		reactantIndex = reactant->getId() - 1;

		// Get the partial derivatives
		reactant->getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];

			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	// Update the column in the Jacobian that represents the moment for the super clusters
	for (int i = 0; i < superClusters.size(); i++) {
		auto reactant = std::static_pointer_cast<PSISuperCluster>(superClusters[i]);

		// Get the super cluster index
		reactantIndex = reactant->getId() - 1;

		// Get the partial derivatives
		reactant->getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the helium momentum index
		reactantIndex = reactant->getHeMomentumId() - 1;

		// Get the partial derivatives
		reactant->getHeMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the vacancy momentum index
		reactantIndex = reactant->getVMomentumId() - 1;

		// Get the partial derivatives
		reactant->getVMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	return;
}

double PSIClusterReactionNetwork::computeBindingEnergy(const DissociationReaction& reaction) const {

	double bindingEnergy = 5.0;
	if (reaction.dissociating->getType() == ReactantType::He
			&& reaction.first->getType() == ReactantType::He) {
		if (reaction.dissociating->getSize() == 2)
			bindingEnergy = 0.5;
		else
			bindingEnergy = 1.0;
	}
	if (reaction.dissociating->getType() == ReactantType::V
			&& reaction.first->getType() == ReactantType::V) {
		int size = reaction.dissociating->getSize();
		bindingEnergy = 1.73
				- 2.59
						* (pow((double) size, 2.0 / 3.0)
								- pow((double) size - 1.0, 2.0 / 3.0));
	}
	if ((reaction.dissociating->getType() == ReactantType::HeV)
			&& (reaction.first->getType() == ReactantType::V
					|| reaction.second->getType() == ReactantType::V)) {
		auto& comp = reaction.dissociating->getComposition();
		bindingEnergy = 1.73
				- 2.59
						* (pow((double) comp[toCompIdx(Species::V)] , 2.0 / 3.0)
								- pow((double) comp[toCompIdx(Species::V)] - 1.0, 2.0 / 3.0))
				+ 2.5
						* log(
								1.0
										+ ((double) comp[toCompIdx(Species::He)]
												/ (double) comp[toCompIdx(Species::V)] ));
	}
	if (reaction.dissociating->getType() == ReactantType::PSISuper
			&& (reaction.first->getType() == ReactantType::V
					|| reaction.second->getType() == ReactantType::V)) {
		auto superCluster = (PSISuperCluster *) reaction.dissociating;
		auto& comp = reaction.dissociating->getComposition();
		double numV = (double) comp[toCompIdx(Species::V)] ;
		double numHe = (double) comp[toCompIdx(Species::He)] ;
		bindingEnergy = 1.73
				- 2.59 * (pow(numV, 2.0 / 3.0) - pow(numV - 1.0, 2.0 / 3.0))
				+ 2.5 * log(1.0 + (numHe / numV));
	}
	if (reaction.first->getType() == ReactantType::I
			|| reaction.second->getType() == ReactantType::I) {
		if (reaction.dissociating->getType() == ReactantType::HeV) {
			auto& comp = reaction.dissociating->getComposition();
			bindingEnergy =
					4.88
							+ 2.59
									* (pow((double) comp[toCompIdx(Species::V)] , 2.0 / 3.0)
											- pow((double) comp[toCompIdx(Species::V)] - 1.0,
													2.0 / 3.0))
							- 2.5
									* log(
											1.0
													+ ((double) comp[toCompIdx(Species::He)]
															/ (double) comp[toCompIdx(Species::V)] ));
		} else if (reaction.dissociating->getType() == ReactantType::PSISuper) {
			auto superCluster = (PSISuperCluster *) reaction.dissociating;
			auto& comp = reaction.dissociating->getComposition();
			double numV = (double) comp[toCompIdx(Species::V)] ;
			double numHe = (double) comp[toCompIdx(Species::He)] ;
			bindingEnergy = 4.88
					+ 2.59 * (pow(numV, 2.0 / 3.0) - pow(numV - 1.0, 2.0 / 3.0))
					- 2.5 * log(1.0 + (numHe / numV));
		} else if (reaction.dissociating->getType() == ReactantType::He) {
			int size = reaction.dissociating->getSize();
			switch (size) {
			case 1:
				bindingEnergy = 4.31;
				break;
			case 2:
				bindingEnergy = 2.90;
				break;
			case 3:
				bindingEnergy = 2.02;
				break;
			case 4:
				bindingEnergy = 1.09;
				break;
			case 5:
				bindingEnergy = 0.58;
				break;
			case 6:
				bindingEnergy = 0.13;
				break;
			case 7:
				bindingEnergy = -0.25;
				break;
			case 8:
				bindingEnergy = -0.59;
				break;
			default:
				break;
			}
		}

	}

//	if (bindingEnergy < -5.0)
//	std::cout << "dissociation: " << reaction.dissociating->getName() << " -> "
//			<< reaction.first->getName() << " + "
//			<< reaction.second->getName() << " : " << bindingEnergy
//			<< std::endl;

	return max(bindingEnergy, -5.0);
}


IReactant * PSIClusterReactionNetwork::getSuperFromComp(IReactant::SizeType nHe, IReactant::SizeType nV) {

    IReactant* ret = nullptr;

    auto heIntervalBase = findBoundsIntervalBase(nHe);
    auto vIntervalBase = findBoundsIntervalBase(nV);

    if((heIntervalBase != 0) and (vIntervalBase != 0)) {
        auto iter = superClusterLookupMap.find(std::make_pair(heIntervalBase, vIntervalBase));
        if(iter != superClusterLookupMap.end()) {
            IReactant& super = iter->second;    // Get ref from reference_wrapper.
            ret = &super;                      // Return pointer to actual object.
            assert(static_cast<PSISuperCluster*>(ret)->isIn(nHe, nV));
        }
    }

    return ret;
}

} // namespace xolotlCore

