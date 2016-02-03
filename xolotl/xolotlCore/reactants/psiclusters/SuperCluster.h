#ifndef SUPERCLUSTER_H
#define SUPERCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>

namespace xolotlCore {

/**
 *  A cluster gathering the average properties of many HeV clusters.
 */
class SuperCluster: public PSICluster {

private:

	//! The mean number of helium atoms in this cluster.
	double numHe;

	//! The mean number of atomic vacancies in this cluster.
	double numV;

	//! The number of HeV clusters gathered in this one.
	int sectionWidth;

	//! The 0th order momentum (mean).
	double l0;

	//! The first order momentum.
	double l1;

	//! The dispersion in the group
	double dispersion;

	//! The map containing all the reacting pairs separated by original helium cluster size.
	std::map <int, std::vector<ClusterPair> > reactingMap;

	//! The map containing all the combining clusters separated by original helium cluster size.
	std::map <int, std::vector<CombiningCluster> > combiningMap;

	//! The map containing all the dissociating pairs separated by original helium cluster size.
	std::map <int, std::vector<ClusterPair> > dissociatingMap;

	//! The map containing all the emission pairs separated by original helium cluster size.
	std::map <int, std::vector<ClusterPair> > emissionMap;

	//! The map containing all the effective reacting pairs separated by original helium cluster size.
	std::map <int, std::vector<ClusterPair *> > effReactingMap;

	//! The map containing all the effective combining clusters separated by original helium cluster size.
	std::map <int, std::vector<CombiningCluster *> > effCombiningMap;

	//! The map containing all the effective dissociating pairs separated by original helium cluster size.
	std::map <int, std::vector<ClusterPair *> > effDissociatingMap;

	//! The map containing all the effective emission pairs separated by original helium cluster size.
	std::map <int, std::vector<ClusterPair *> > effEmissionMap;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	SuperCluster() :
		PSICluster(1)
	{ numHe = 1.0; numV = 1.0; }

public:

	//! The vector of HeV clusters it will replace
	std::vector<PSICluster *> heVVector;

	/**
	 * The constructor. All SuperClusters must be initialized with its
	 * composition.
	 *
	 * @param numHe The mean number of helium atoms in this cluster
	 * @param numV The mean number of vacancies in this cluster
	 * @param width The width of this super cluster
	 * @param radius The mean radius
	 * @param energy The mean formation energy
	 * @param registry The performance handler registry
	 */
	SuperCluster(double numHe, double numV, int width, double radius, double energy,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor.
	 *
	 * @param other the reactant to be copied
	 */
	SuperCluster(const SuperCluster &other);

	//! Destructor
	~SuperCluster() {}

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of SuperCluster.
	 *
	 * @return A copy of this reactant
	 */
	virtual std::shared_ptr<Reactant> clone();

	/**
	 * This operation returns true to signify that this cluster is a mixture of
	 * He and V.
	 *
	 * @return True if mixed
	 */
	virtual bool isMixed() const {return true;}

	/**
	 * Set the HeV vector
	 */
	void setHeVVector(std::vector<PSICluster *> vec)
		{heVVector = vec;}

	/**
	 * This operation returns the current concentration.
	 *
	 * @param id The id in the group
	 * @return The concentration of this reactant
	 */
	double getConcentration(double id) const;

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of helium
	 * @return The distance to the mean number of helium in the group
	 */
	double getDistance(int he) const;

	/**
	 * Computes a row of the reaction connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants alone can form a reaction, the element at the position
	 * of the second reactant is 1, otherwise 0.
	 */
	void createReactionConnectivity();

	/**
	 * Computes a row of the dissociation connectivity matrix corresponding to
	 * this reactant.
	 *
	 * If two reactants together can be produced by a single reaction,
	 * the element at the position of the second reactant is 1, otherwise 0.
	 */
	void createDissociationConnectivity();

	/**
	 * Calculate all the rate constants for the reactions and dissociations in which this
	 * cluster is taking part. Store these values in the kConstant field of ClusterPair
	 * or CombiningCluster. Need to be called only when the temperature changes.
	 */
	void computeRateConstants();

	/**
	 * This operation sets the zeroth order momentum.
	 *
	 * @param mom The momentum
	 */
	void setZerothMomentum(double mom) {l0 = mom;}

	/**
	 * This operation sets the first order momentum.
	 *
	 * @param mom The momentum
	 */
	void setFirstMomentum(double mom) {l1 = mom;}

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the production and dissociation vectors.
	 */
	void resetConnectivities();

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it.
	 *
	 * @return The flux due to dissociation of other clusters
	 */
	double getDissociationFlux() const;

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux() const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters.
	 *
	 * @return The flux due to this cluster being produced
	 */
	double getProductionFlux() const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others.
	 *
	 * @return The flux due to this cluster combining with other clusters
	 */
	double getCombinationFlux() const;

	/**
	 * This operation returns the moment flux of this cluster in the
	 * current network.
	 *
	 * @return The total change in the moment flux
	 */
	double getMomentFlux() const;

	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getProductionPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getCombinationPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getDissociationPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 */
	void getEmissionPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives for the momentum.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getMomentPartialDerivatives(std::vector<double> & partials) const;

};
//end class SuperCluster

} /* end namespace xolotlCore */
#endif
