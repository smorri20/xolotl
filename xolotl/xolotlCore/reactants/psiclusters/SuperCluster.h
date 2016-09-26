#ifndef SUPERCLUSTER_H
#define SUPERCLUSTER_H

// Includes
#include "PSICluster.h"
#include <string>
#include <forward_list>

namespace xolotlCore {

/**
 *  A cluster gathering the average properties of many HeV clusters.
 */
class SuperCluster: public PSICluster {

protected:

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two body production reactions.
	 *
	 * The constants are stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	class SuperClusterProductionPair {
	public:

		/**
		 * The first cluster in the pair
		 */
		PSICluster * first;

		/**
		 * The second cluster in the pair
		 */
		PSICluster * second;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		double kConstant;

		/**
		 * All the coefficient needed to compute each element
		 */
		double a00000;
		double a00001;
		double a00010;
		double a00011;
		double a00100;
		double a00101;
		double a00110;
		double a00111;
		double a01000;
		double a01001;
		double a01010;
		double a01011;
		double a01100;
		double a01101;
		double a01110;
		double a01111;
		double a10000;
		double a10001;
		double a10010;
		double a10011;
		double a10100;
		double a10101;
		double a10110;
		double a10111;
		double a11000;
		double a11001;
		double a11010;

		//! The constructor
		SuperClusterProductionPair(PSICluster * firstPtr, PSICluster * secondPtr, double k)
		: first(firstPtr), second(secondPtr), kConstant(k), a00000(0.0), a00001(0.0), a00010(0.0),
		  a00011(0.0), a00100(0.0), a00101(0.0), a00110(0.0), a00111(0.0), a01000(0.0), a01001(0.0),
		  a01010(0.0), a01011(0.0), a01100(0.0), a01101(0.0), a01110(0.0), a01111(0.0), a10000(0.0),
		  a10001(0.0), a10010(0.0), a10011(0.0), a10100(0.0), a10101(0.0), a10110(0.0), a10111(0.0),
		  a11000(0.0), a11001(0.0), a11010(0.0) {}
	};

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two dissociation reactions.
	 *
	 * The constants are stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	class SuperClusterDissociationPair {
	public:

		/**
		 * The first cluster in the pair
		 */
		PSICluster * first;

		/**
		 * The second cluster in the pair
		 */
		PSICluster * second;

		/**
		 * The reaction/dissociation constant associated to this
		 * reaction or dissociation
		 */
		double kConstant;

		/**
		 * All the coefficient needed to compute each element
		 */
		double a0000;
		double a0001;
		double a0010;
		double a0011;
		double a0100;
		double a0101;
		double a0110;
		double a0111;
		double a1000;

		//! The constructor
		SuperClusterDissociationPair(PSICluster * firstPtr, PSICluster * secondPtr, double k)
		: first(firstPtr), second(secondPtr), kConstant(k), a0000(0.0), a0001(0.0), a0010(0.0),
		  a0011(0.0), a0100(0.0), a0101(0.0), a0110(0.0), a0111(0.0), a1000(0.0) {}
	};

private:

	//! The mean number of helium atoms in this cluster.
	double numHe;

	//! The mean number of atomic vacancies in this cluster.
	double numV;

	//! The total number of clusters gathered in this super cluster.
	int nTot;

	//! The width in the helium direction.
	int sectionHeWidth;

	//! The width in the vacancy direction.
	int sectionVWidth;

	//! The 0th order momentum (mean).
	double l0;

	//! The first order momentum in the helium direction.
	double l1He;

	//! The first order momentum in the vacancy direction.
	double l1V;

	//! The dispersion in the group in the helium direction.
	double dispersionHe;

	//! The dispersion in the group in the vacancy direction.
	double dispersionV;

	//! The map containing all the reacting pairs separated by original composition.
	std::map <std::pair<int, int>, std::vector<ClusterPair> > reactingMap;

	//! The map containing all the combining clusters separated by original composition.
	std::map <std::pair<int, int>, std::vector<CombiningCluster> > combiningMap;

	//! The map containing all the dissociating pairs separated by original composition.
	std::map <std::pair<int, int>, std::vector<ClusterPair> > dissociatingMap;

	//! The map containing all the emission pairs separated by original composition.
	std::map <std::pair<int, int>, std::vector<ClusterPair> > emissionMap;

	//! The map containing all the effective reacting pairs separated by original composition.
	std::map <std::pair<int, int>, std::vector<ClusterPair *> > effReactingMap;

	//! The map containing all the effective combining clusters separated by original composition.
	std::map <std::pair<int, int>, std::vector<CombiningCluster *> > effCombiningMap;

	//! The map containing all the effective dissociating pairs separated by original composition.
	std::map <std::pair<int, int>, std::vector<ClusterPair *> > effDissociatingMap;

	//! The map containing all the effective emission pairs separated by original composition.
	std::map <std::pair<int, int>, std::vector<ClusterPair *> > effEmissionMap;

	//! The list of optimized effective reacting pairs.
	std::forward_list <SuperClusterProductionPair> effReactingList;

	//! The list of optimized effective combining pairs.
	std::forward_list <SuperClusterProductionPair> effCombiningList;

	//! The list of optimized effective dissociating pairs.
	std::forward_list <SuperClusterDissociationPair> effDissociatingList;

	//! The list of optimized effective emission pairs.
	std::forward_list <SuperClusterDissociationPair> effEmissionList;

	/**
	 * The vector containing a boolean to know if this super cluster burst at each depth.
	 */
	std::vector<bool> burstingIndexVector;

	/**
	 * The default constructor is private because PSIClusters must always be
	 * initialized with a size.
	 */
	SuperCluster() :
		PSICluster() {}

	/**
	 * Group the same reactions together.
	 */
	void optimizeReactions();

public:

	//! The vector of HeV clusters it will replace
	std::vector<PSICluster *> heVVector;

	/**
	 * The constructor. All SuperClusters must be initialized with its
	 * composition.
	 *
	 * @param numHe The mean number of helium atoms in this cluster
	 * @param numV The mean number of vacancies in this cluster
	 * @param nTot The total number of clusters in this cluster
	 * @param heWidth The width of this super cluster in the helium direction
	 * @param vWidth The width of this super cluster in the vacancy direction
	 * @param radius The mean radius
	 * @param registry The performance handler registry
	 */
	SuperCluster(double numHe, double numV, int nTot, int heWidth, int vWidth,  double radius,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * Copy constructor.
	 *
	 * @param other the reactant to be copied
	 */
	SuperCluster(SuperCluster &other);

	//! Destructor
	~SuperCluster() {}

	/**
	 * This operation returns a Reactant that is created using the copy
	 * constructor of SuperCluster.
	 *
	 * @return A copy of this reactant
	 */
	virtual std::shared_ptr<IReactant> clone();

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
	 * @param distHe The helium distance in the group
	 * @param distV The vacancy distance in the group
	 * @return The concentration of this reactant
	 */
	double getConcentration(double distHe, double distV) const;

	/**
	 * This operation returns the first helium momentum.
	 *
	 * @return The momentum
	 */
	double getHeMomentum() const;

	/**
	 * This operation returns the first vacancy momentum.
	 *
	 * @return The momentum
	 */
	double getVMomentum() const;

	/**
	 * This operation returns the current total concentration of clusters in the group.

	 * @return The concentration
	 */
	double getTotalConcentration() const;

	/**
	 * This operation returns the current total concentration of helium in the group.

	 * @return The concentration
	 */
	double getTotalHeliumConcentration() const;

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of helium
	 * @return The distance to the mean number of helium in the group
	 */
	double getHeDistance(int he) const;

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param he The number of vacancy
	 * @return The distance to the mean number of vacancy in the group
	 */
	double getVDistance(int v) const;

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
	 * or CombiningCluster. Need to be called only the first time.
	 */
	void computeRateConstants();

	/**
	 * Calculate all the rate constants for the reactions and dissociations in which this
	 * cluster is taking part. Store these values in the kConstant field of ClusterPair
	 * or CombiningCluster. Need to be called only when the temperature changes.
	 */
	void updateRateConstants();

	/**
	 * This operation sets the zeroth order momentum.
	 *
	 * @param mom The momentum
	 */
	void setZerothMomentum(double mom) {l0 = mom;}

	/**
	 * This operation sets the first order momentum in the helium direction.
	 *
	 * @param mom The momentum
	 */
	void setHeMomentum(double mom) {l1He = mom;}

	/**
	 * This operation sets the first order momentum in the vacancy direction.
	 *
	 * @param mom The momentum
	 */
	void setVMomentum(double mom) {l1V = mom;}

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the production and dissociation vectors.
	 */
	void resetConnectivities();

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to dissociation of other clusters
	 */
	double getDissociationFlux();

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to its dissociation
	 */
	double getEmissionFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to this cluster being produced
	 */
	double getProductionFlux();

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others. Compute the contributions to
	 * the momentum fluxes at the same time.
	 *
	 * @return The flux due to this cluster combining with other clusters
	 */
	double getCombinationFlux();

	/**
	 * This operation works as getPartialDerivatives above, but instead of
	 * returning a vector that it creates it fills a vector that is passed to
	 * it by the caller. This allows the caller to optimize the amount of
	 * memory allocations to just one if they are accessing the partial
	 * derivatives many times.
	 *
	 * @param the vector that should be filled with the partial derivatives
	 * for this reactant where index zero corresponds to the first reactant in
	 * the list returned by the ReactionNetwork::getAll() operation. The size of
	 * the vector should be equal to ReactionNetwork::size().
	 *
	 */
	virtual void getPartialDerivatives(std::vector<double> & partials) const;

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
	 * This operation computes the partial derivatives for the helium momentum.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getHeMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives for the vacancy momentum.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted.
	 */
	void getVMomentPartialDerivatives(std::vector<double> & partials) const;

	/**
	 * Returns the average number of vacancies.
	 *
	 * @return The average number of vacancies
	 */
	double getNumV() {return numV;}

};
//end class SuperCluster

} /* end namespace xolotlCore */
#endif
