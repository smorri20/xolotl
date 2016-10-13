#include "NEClusterReactionNetwork.h"
#include "NECluster.h"
#include <xolotlPerf.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Constants.h>
#include <NESuperCluster.h>

using namespace xolotlCore;

void NEClusterReactionNetwork::setDefaultPropsAndNames() {
	// Shared pointers for the cluster type map
	std::shared_ptr<std::vector<std::shared_ptr<IReactant>>>xeVector =
	std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> vVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> iVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> xeVVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> xeIVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();
	std::shared_ptr < std::vector<std::shared_ptr<IReactant>>> superVector
	= std::make_shared<std::vector<std::shared_ptr<IReactant>>>();

	// Initialize default properties
	(*properties)["reactionsEnabled"] = "true";
	(*properties)["dissociationsEnabled"] = "true";
	(*properties)["numXeClusters"] = "0";
	(*properties)["numVClusters"] = "0";
	(*properties)["numIClusters"] = "0";
	(*properties)["numXeVClusters"] = "0";
	(*properties)["numXeIClusters"] = "0";
	(*properties)["numSuperClusters"] = "0";
	(*properties)["maxXeClusterSize"] = "0";
	(*properties)["maxVClusterSize"] = "0";
	(*properties)["maxIClusterSize"] = "0";
	(*properties)["maxXeVClusterSize"] = "0";
	(*properties)["maxXeIClusterSize"] = "0";

	// Initialize the current and last size to 0
	networkSize = 0;
	// Set the reactant names
	names.push_back(xeType);
	names.push_back(vType);
	names.push_back(iType);
	names.push_back(NESuperType);
	// Set the compound reactant names
	compoundNames.push_back(xeVType);
	compoundNames.push_back(xeIType);

	// Setup the cluster type map
	clusterTypeMap[xeType] = xeVector;
	clusterTypeMap[vType] = vVector;
	clusterTypeMap[iType] = iVector;
	clusterTypeMap[xeVType] = xeVVector;
	clusterTypeMap[xeIType] = xeIVector;
	clusterTypeMap[NESuperType] = superVector;

	return;
}

NEClusterReactionNetwork::NEClusterReactionNetwork() :
		ReactionNetwork() {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

NEClusterReactionNetwork::NEClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork(registry) {
	// Setup the properties map and the name lists
	setDefaultPropsAndNames();

	return;
}

NEClusterReactionNetwork::NEClusterReactionNetwork(
		const NEClusterReactionNetwork &other) :
		ReactionNetwork(other) {
	// The size and ids do not need to be copied. They will be fixed when the
	// reactants are added.

	// Reset the properties table so that it can be properly updated when the
	// network is filled.
	setDefaultPropsAndNames();
	// Get all of the reactants from the other network and add them to this one
	// Load the single-species clusters. Calling getAll() will not work because
	// it is not const.
	std::vector<std::shared_ptr<IReactant> > reactants;
	for (auto it = other.singleSpeciesMap.begin();
			it != other.singleSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the mixed-species clusters
	for (auto it = other.mixedSpeciesMap.begin();
			it != other.mixedSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	// Load the super-species clusters
	for (auto it = other.superSpeciesMap.begin();
			it != other.superSpeciesMap.end(); ++it) {
		reactants.push_back(it->second);
	}
	for (unsigned int i = 0; i < reactants.size(); i++) {
		add(reactants[i]->clone());
	}

	return;
}

void NEClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	for (int i = 0; i < networkSize; i++) {
		// Now that the diffusion coefficients of all the reactants
		// are updated, the reaction and dissociation rates can be
		// recomputed
		auto cluster = allReactants->at(i);
		cluster->computeRateConstants();
	}

	return;
}

double NEClusterReactionNetwork::getTemperature() const {
	return temperature;
}

IReactant * NEClusterReactionNetwork::get(const std::string& type,
		const int size) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { xeType, 0 },
			{ vType, 0 }, { iType, 0 }, { heType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Setup the composition map to default values because it is static
	composition[xeType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name and size are valid
	if ((type == xeType || type == vType || type == iType) && size >= 1) {
		composition[type] = size;
		//std::string encodedName = NECluster::encodeCompositionAsName(composition);
		// Make sure the reactant is in the map
        std::string compstr = Reactant::toCanonicalString(composition);
		if (singleSpeciesMap.count(compstr)) {
			retReactant = singleSpeciesMap.at(compstr);
		}
	}

	return retReactant.get();
}

IReactant * NEClusterReactionNetwork::getCompound(const std::string& type,
		const std::vector<int>& sizes) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { xeType, 0 },
			{ vType, 0 }, { iType, 0 }, { heType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Setup the composition map to default values because it is static
	composition[xeType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if ((type == xeVType || type == xeIType) && sizes.size() == 3) {
		composition[xeType] = sizes[0];
		composition[vType] = sizes[1];
		composition[iType] = sizes[2];
		// Make sure the reactant is in the map
        std::string compstr = Reactant::toCanonicalString(composition);
		if (mixedSpeciesMap.count(compstr)) {
			retReactant = mixedSpeciesMap.at(compstr);
		}
	}

	return retReactant.get();
}

IReactant * NEClusterReactionNetwork::getSuper(const std::string& type,
		const int size) const {
	// Local Declarations
	static std::map<std::string, int> composition = { { xeType, 0 },
			{ vType, 0 }, { iType, 0 }, { heType, 0 } };
	std::shared_ptr<IReactant> retReactant;

	// Setup the composition map to default values
	composition[xeType] = 0;
	composition[vType] = 0;
	composition[iType] = 0;

	// Only pull the reactant if the name and size are valid.
	if (type == NESuperType && size >= 1) {
		composition[xeType] = size;
		// Make sure the reactant is in the map
        std::string compstr = Reactant::toCanonicalString(composition);
		if (superSpeciesMap.count(compstr)) {
			retReactant = superSpeciesMap.at(compstr);
		}
	}

	return retReactant.get();
}

const std::shared_ptr<std::vector<IReactant *>> & NEClusterReactionNetwork::getAll() const {
	return allReactants;
}

std::vector<IReactant *> NEClusterReactionNetwork::getAll(
		const std::string& name) const {
	// Local Declarations
	std::vector<IReactant *> reactants;

	// Only pull the reactants if the name is valid
	if (name == xeType || name == vType || name == iType || name == xeVType
			|| name == xeIType || name == NESuperType) {
		std::shared_ptr<std::vector<std::shared_ptr<IReactant>> > storedReactants =
				clusterTypeMap.at(name);
		int vecSize = storedReactants->size();
		for (int i = 0; i < vecSize; i++) {
			reactants.push_back(storedReactants->at(i).get());
		}
	}

	return reactants;
}

void NEClusterReactionNetwork::add(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numXe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	std::string numClusterKey, clusterSizeKey;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto composition = reactant->getComposition();
        std::string compstr = reactant->getCompositionString();

		// Get the species sizes
		numXe = composition.at(xeType);
		numV = composition.at(vType);
		numI = composition.at(iType);

		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numXe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
		if (isMixed && mixedSpeciesMap.count(compstr) == 0) {
			// Put the compound in its map
			mixedSpeciesMap[compstr] = reactant;
			// Figure out whether we have XeV or XeI and set the keys
			if (numV > 0) {
				numClusterKey = "numXeVClusters";
				clusterSizeKey = "maxXeVClusterSize";
			} else {
				numClusterKey = "numXeIClusters";
				clusterSizeKey = "maxXeIClusterSize";
			}
		} else if (!isMixed && singleSpeciesMap.count(compstr) == 0) {
			/// Put the reactant in its map
			singleSpeciesMap[compstr] = reactant;

			// Figure out whether we have Xe, V or I and set the keys
			if (numXe > 0) {
				numClusterKey = "numXeClusters";
				clusterSizeKey = "maxXeClusterSize";
			} else if (numV > 0) {
				numClusterKey = "numVClusters";
				clusterSizeKey = "maxVClusterSize";
			} else {
				numClusterKey = "numIClusters";
				clusterSizeKey = "maxIClusterSize";
			}
		} else {
			std::stringstream errStream;
			errStream << "NEClusterReactionNetwork Message: "
					<< "Duplicate Reactant (Xe=" << numXe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the number of total clusters of this type
		int numClusters = std::stoi(properties->at(numClusterKey));
		numClusters++;
		(*properties)[numClusterKey] = std::to_string((long long) numClusters);
		// Increment the max cluster size key
		int maxSize = std::stoi(properties->at(clusterSizeKey));
		int clusterSize = numXe + numV + numI;
		maxSize = std::max(clusterSize, maxSize);
		(*properties)[clusterSizeKey] = std::to_string((long long) maxSize);
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);
		// Get the vector for this reactant from the type map
		auto clusters = clusterTypeMap[reactant->getType()];

		clusters->push_back(reactant);
		// Add the pointer to the list of all clusters
		allReactants->push_back(reactant.get());
	}

	return;
}

void NEClusterReactionNetwork::addSuper(std::shared_ptr<IReactant> reactant) {
	// Local Declarations
	int numXe = 0, numV = 0, numI = 0;
	bool isMixed = false;
	std::string numClusterKey;

	// Only add a complete reactant
	if (reactant != NULL) {
		// Get the composition
		auto composition = reactant->getComposition();
		// Get the species sizes
		numXe = composition.at(xeType);
		numV = composition.at(vType);
		numI = composition.at(iType);
		// Determine if the cluster is a compound. If there is more than one
		// type, then the check below will sum to greater than one and we know
		// that we have a mixed cluster.
		isMixed = ((numXe > 0) + (numV > 0) + (numI > 0)) > 1;
		// Only add the element if we don't already have it
		// Add the compound or regular reactant.
        std::string compstr = reactant->getCompositionString();
		if (!isMixed && superSpeciesMap.count(compstr) == 0) {
			// Put the compound in its map
			superSpeciesMap[compstr] = reactant;
			// Set the key
			numClusterKey = "numSuperClusters";
		} else {
			std::stringstream errStream;
			errStream << "NEClusterReactionNetwork Message: "
					<< "Duplicate Super Reactant (Xe=" << numXe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}

		// Increment the number of total clusters of this type
		int numClusters = std::stoi(properties->at(numClusterKey));
		numClusters++;
		(*properties)[numClusterKey] = std::to_string((long long) numClusters);
		// Update the size
		++networkSize;
		// Set the id for this cluster
		reactant->setId(networkSize);
		// Get the vector for this reactant from the type map
		auto clusters = clusterTypeMap[reactant->getType()];
		clusters->push_back(reactant);
		// Add the pointer to the list of all clusters
		allReactants->push_back(reactant.get());
	}

	return;
}

void NEClusterReactionNetwork::removeReactant(IReactant * reactant) {
	auto comp = reactant->getComposition();
	std::string type = reactant->getType();

	// Look for the reactant in allReactants
	for (auto it = allReactants->begin(); it != allReactants->end(); ++it) {
		auto tempType = (*it)->getType();
		// Compare the types and skip if necessary
		if (tempType != type) continue;

		auto tempComp = (*it)->getComposition();
		// Compare the compositions and skip if necessary
		if (tempComp[xeType] != comp[xeType] || tempComp[vType] != comp[vType]
						|| tempComp[iType] != comp[iType]) continue;

		allReactants->erase(it);
		break;
	}

	// Look for the reactant in clusterTypeMap
	auto clusters = clusterTypeMap[type];
	for (auto it = clusters->begin(); it != clusters->end(); ++it) {
		auto tempType = (*it)->getType();
		// Compare the types and skip if necessary
		if (tempType != type) continue;

		auto tempComp = (*it)->getComposition();
		// Compare the compositions and skip if necessary
		if (tempComp[xeType] != comp[xeType] || tempComp[vType] != comp[vType]
						|| tempComp[iType] != comp[iType]) continue;

		clusters->erase(it);
		break;
	}

	// Look for the reactant in mixedSpeciesMap
	mixedSpeciesMap.erase(reactant->getCompositionString());

	return;
}

void NEClusterReactionNetwork::reinitializeNetwork() {
	// Reset the Ids
	int id = 0;
	for (auto it = allReactants->begin(); it != allReactants->end(); ++it) {
		id++;
		(*it)->setId(id);
		(*it)->setMomentumId(id);
	}

	// Reset the network size
	networkSize = id;

	// Get all the super clusters and loop on them
	for (auto it = clusterTypeMap[NESuperType]->begin(); it != clusterTypeMap[NESuperType]->end(); ++it) {
		id++;
		(*it)->setMomentumId(id);
	}

	return;
}

void NEClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
	for (auto it = allReactants->begin(); it != allReactants->end(); ++it) {
		(*it)->resetConnectivities();
	}

	return;
}

void NEClusterReactionNetwork::setProperty(const std::string& key,
		const std::string& value) {
	// Check the keys and value before trying to set the property
	if (!key.empty() && !value.empty() && key != "numXeClusters"
			&& key != "numVClusters" && key != "numIClusters"
			&& key != "numSuperClusters" && key != "maxXeClusterSize"
			&& key != "maxVClusterSize" && key != "maxIClusterSize"
			&& key != "maxXeVClusterSize" && key != "maxXeIClusterSize") {
		// Add the property if it made it through that!
		(*properties)[key] = value;
	}

	return;
}

void NEClusterReactionNetwork::updateConcentrationsFromArray(double * concentrations) {
	// Local Declarations
	auto reactants = getAll();
	int size = reactants->size();
	int id = 0;

	// Set the concentrations
	concUpdateCounter->increment();	// increment the update concentration counter
	for (int i = 0; i < size; i++) {
		id = reactants->at(i)->getId() - 1;
		reactants->at(i)->setConcentration(concentrations[id]);
	}

	// Set the moments
	int numSuperClusters = stoi(properties->at("numSuperClusters"));
	for (int i = size - numSuperClusters; i < size; i++) {
		// Get the superCluster
		auto cluster = (NESuperCluster *) reactants->at(i);
		id = cluster->getId() - 1;
		cluster->setZerothMomentum(concentrations[id]);
		id = cluster->getMomentumId() - 1;
		cluster->setMomentum(concentrations[id]);
	}

	return;
}


