#ifndef PSICLUSTER_H
#define PSICLUSTER_H

// Includes
#include <Reactant.h>
#include "IntegerRange.h"
#include "NDArray.h"

namespace xolotlPerf {
class ITimer;
}

namespace xolotlCore {

class PSIClusterReactionNetwork;

/**
 * The PSICluster class is a Reactant that is specialized to work for
 * simulations of plasma-surface interactions. It provides special routines
 * for calculating the total flux due to production and dissociation and
 * obtaining the cluster size.
 *
 * PSIClusters must always be initialized with a size. If the constructor is
 * passed a size of zero or less, the actual size will be set to 1.
 *
 * The getComposition() operation is implemented by subclasses and will always
 * return a map with the keys He, V, I, HeV or HeI. The operation getTypeName()
 * will always return one of the same values.
 *
 * As a rule, it is possible to access directly some of the private members of
 * this class (id, concentration, reactionRadius, diffusionCoefficient, size,
 * type) instead of using the "get" functions for performance reasons. In
 * order to change these values the "set" functions must still be used.
 */
class PSICluster: public Reactant {

public:
    // Concise name for type used for the PSI index list.
    using IndexList = Array<int, 5>;

protected:

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for two body reactions or dissociation.
	 *
	 * The constant k+ or k- is stored along the clusters taking part in the
	 * reaction or dissociation for faster computation because they only change
	 * when the temperature change. k is computed when setTemperature() is called.
	 */
	struct ClusterPairBase {

		/**
		 * The first cluster in the pair
		 */
		PSICluster& first;

		/**
		 * The second cluster in the pair
		 */
		PSICluster& second;

		/**
		 * The reaction/dissociation pointer to the list
		 */
		Reaction& reaction;

		//! The constructor
		ClusterPairBase(Reaction& _reaction,
                PSICluster& _first, PSICluster& _second) :
				first(_first), second(_second), reaction(_reaction) {
		}

		/**
		 * Default constructor, disallowed.
		 */
		ClusterPairBase() = delete;

		// NB: if PSICluster keeps these in a std::vector,
		// copy ctor is needed.
		ClusterPairBase(const ClusterPairBase& other) :
                first(other.first), second(other.second),
                reaction(other.reaction) {
		}
	};

	struct ClusterPair : public ClusterPairBase {

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the moment of A, the second of B
		 * in A + B -> C
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> D
		 * 3 -> T
		 * 4 -> V
		 */
        Array<double, 5, 5> coefs;

		//! The constructor
		ClusterPair(Reaction& _reaction,
                PSICluster& _first, PSICluster& _second) :
            ClusterPairBase(_reaction, _first, _second) {

            coefs.Init(0);
		}

		/**
		 * Default constructor, disallowed.
		 */
		ClusterPair() = delete;

		// NB: if PSICluster keeps these in a std::vector,
		// copy ctor is needed.
		ClusterPair(const ClusterPair& other) :
            ClusterPairBase(other),
            coefs(other.coefs) {
		}
	};

	struct ClusterPair0 : public ClusterPairBase {

        /**
         * Scalar version of coefs[0][0] supporting
         * faster reads for zeroth-moment-only flux/partials computations.
         */
        double coeff0;

		//! The constructor
		ClusterPair0(Reaction& _reaction,
                PSICluster& _first, PSICluster& _second) :
            ClusterPairBase(_reaction, _first, _second),
            coeff0(0) {
		}

		/**
		 * Default constructor, disallowed.
		 */
		ClusterPair0() = delete;

        /**
         * Construct as copy of given cluster.
         */
		ClusterPair0(const ClusterPair0& other) :
            ClusterPairBase(other),
            coeff0(other.coeff0) {
		}

        /**
         * Construct as copy of given cluster with full set of coefficients.
         */
        ClusterPair0(const ClusterPair& other) :
            ClusterPairBase(other),
            coeff0(other.coefs[0][0]) {
        }
	};

	/**
	 * This is a protected class that is used to implement the flux calculations
	 * for combinations.
	 *
	 * The constant k+ is stored along the cluster that combines with this cluster
	 * for faster computation because they only change when the temperature change.
	 * k+ is computed when setTemperature() is called.
	 */
	struct CombiningClusterBase {

		/**
		 * The combining cluster
		 */
		PSICluster& combining;

		/**
		 * The reaction pointer to the list
		 */
		Reaction& reaction;

		//! The constructor
		CombiningClusterBase(Reaction& _reaction, PSICluster& _comb) :
				combining(_comb), reaction(_reaction) {

		}

		/**
		 * Default constructor, disallowed to prohibit building without args.
		 */
		CombiningClusterBase() = delete;

		// NB: if PSICluster keeps these in a std::vector,
		// copy ctor is needed.
		CombiningClusterBase(const CombiningClusterBase& other) :
                combining(other.combining), reaction(other.reaction) {
		}
	};

	struct CombiningCluster : public CombiningClusterBase {

		/**
		 * All the coefficient needed to compute each element
		 * The first number represent the moment of A
		 * in A + this -> C
		 *
		 * 0 -> l0
		 * 1 -> He
		 * 2 -> D
		 * 3 -> T
		 * 4 -> V
		 */
        Array<double, 5> coefs;

		//! The constructor
		CombiningCluster(Reaction& _reaction, PSICluster& _comb) :
            CombiningClusterBase(_reaction, _comb) {

            coefs.Init(0);
		}

		/**
		 * Default constructor, disallowed to prohibit building without args.
		 */
		CombiningCluster() = delete;

		// NB: if PSICluster keeps these in a std::vector,
		// copy ctor is needed.
		CombiningCluster(const CombiningCluster& other) :
            CombiningClusterBase(other),
			coefs(other.coefs) {
		}
	};

	struct CombiningCluster0 : public CombiningClusterBase {

        /**
         * Scalar version of coefs[0] supporting
         * faster reads for zeroth-moment-only flux/partials computations.
         */
        double coeff0;

		//! The constructor
		CombiningCluster0(Reaction& _reaction, PSICluster& _comb) :
            CombiningClusterBase(_reaction, _comb),
            coeff0(0) {
		}

		/**
		 * Default constructor, disallowed to prohibit building without args.
		 */
		CombiningCluster0() = delete;

		// NB: if PSICluster keeps these in a std::vector,
		// copy ctor is needed.
		CombiningCluster0(const CombiningCluster0& other) :
            CombiningClusterBase(other),
            coeff0(other.coeff0) {
		}

		CombiningCluster0(const CombiningCluster& other) :
            CombiningClusterBase(other),
            coeff0(other.coefs[0]) {
		}
	};

	/**
	 * Bounds on number of He atoms represented by this cluster.
	 */
	IntegerRange<IReactant::SizeType> bounds[4] = { IntegerRange<
			IReactant::SizeType>(0, 0), IntegerRange<IReactant::SizeType>(0, 0),
			IntegerRange<IReactant::SizeType>(0, 0), IntegerRange<
					IReactant::SizeType>(0, 0) };

	//! The dimension of the phase space
	const int& psDim;

	//! The indexList.
    const IndexList& indexList;

	/**
	 * This operation returns a set that contains only the entries of the
	 * reaction connectivity array that are non-zero.
	 *
	 * @return The set of connected reactants. Each entry in the set is the id
	 * of a connected cluster for forward reactions.
	 */
	std::set<int> getReactionConnectivitySet() const;

	/**
	 * This operation returns a set that contains only the entries of the
	 * dissociation connectivity array that are non-zero.
	 *
	 * @return The set of connected reactants. Each entry in the set is the id
	 * of a connected cluster for dissociation reactions
	 */
	const std::set<int> & getDissociationConnectivitySet() const;

	/**
	 * Output coefficients for a given reaction to the given output stream.
	 *
	 * @param os The output stream on which to write the coefficients.
	 * @param curr Information about our participation in a reaction.
	 */
	void dumpCoefficients(std::ostream& os, ClusterPair const& curr) const;
	void dumpCoefficients(std::ostream& os, CombiningCluster const& curr) const;


	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
     * @param[out] flux The flux due to dissociation of other clusters.
	 */
    void getDissociationFlux(const double* __restrict concs, int i,
                                Reactant::Flux& flux) const override;
    void computeDissFlux0(const double* __restrict concs, int xi,
                                Reactant::Flux& flux) const;

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to its dissociation
	 */
    void getEmissionFlux(const double* __restrict concs, int i,
                                Reactant::Flux& flux) const override;
    void computeEmitFlux0(const double* __restrict concs, int xi,
                                Reactant::Flux& flux) const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux flux The flux due to this cluster being produced
	 */
    void getProductionFlux(const double* __restrict concs, int i,
                                Reactant::Flux& flux) const override;
    void computeProdFlux0(const double* __restrict concs,
                                int xi, Reactant::Flux& flux) const;

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to this cluster combining with other clusters
	 */
    void getCombinationFlux(const double* __restrict concs, int i,
                                Reactant::Flux& flux) const override;
    void computeCombFlux0(const double* __restrict concs, int xi,
                                Reactant::Flux& flux) const;


    template<typename FluxType>
    FluxType getTotalFluxHelper0(const double* __restrict concs, int xi) const {

        // Compute the individual fluxes.
        //
        // NOTE: We would much prefer to have the get*Flux() methods 
        // return our type of flux, but our class hierarchy needs them to be 
        // virtual and we have differing Flux types than our base class
        // that we would use a return type.
        FluxType prodFlux;
        computeProdFlux0(concs, xi, prodFlux);

        FluxType combFlux;
        computeCombFlux0(concs, xi, combFlux);

        FluxType dissFlux;
        computeDissFlux0(concs, xi, dissFlux);

        FluxType emitFlux;
        computeEmitFlux0(concs, xi, emitFlux);

        // Compute the total flux.
        return prodFlux - combFlux + dissFlux - emitFlux;
    }

	virtual void computeAllProdPartials0(const double* __restrict concs,
            int xi, std::vector<double> & partials) const;
	void computeOneProdPartials0(const double* __restrict concs,
            int xi, std::vector<double>& partials,
            const ClusterPair0& currPair) const;

	virtual void computeAllCombPartials0(const double* __restrict concs,
            int xi, std::vector<double> & partials) const;
	void computeOneCombPartials0(const double* __restrict concs,
            int xi, std::vector<double>& partials,
            const CombiningCluster0& currCluster) const;

	virtual void computeAllDissPartials0(int xi, 
            std::vector<double> & partials) const;
	void computeOneDissPartials0(int xi, 
            std::vector<double>& partials,
            const ClusterPair0& currPair) const;

	virtual void computeAllEmitPartials0(int xi,
            std::vector<double>& partials) const;



public:

	/**
	 * A vector of ClusterPairs that represents reacting pairs of clusters
	 * that produce this cluster. This vector should be populated early in the
	 * cluster's lifecycle by subclasses. In the standard Xolotl clusters,
	 * this vector is filled in createReactionConnectivity.
	 */
	std::vector<ClusterPair> reactingPairs;
	std::vector<ClusterPair0> reactingPairs0;

	/**
	 * A vector of clusters that combine with this cluster to produce other
	 * clusters. This vector should be populated early in the cluster's
	 * lifecycle by subclasses. In the standard Xolotl clusters, this vector is
	 * filled in createReactionConnectivity.
	 */
	std::vector<CombiningCluster> combiningReactants;
	std::vector<CombiningCluster0> combiningReactants0;

	/**
	 * A vector of pairs of clusters: the first one is the one dissociation into
	 * this cluster, the second one is the one that is emitted at the same time
	 * during the dissociation. This vector should be populated early in the
	 * cluster's lifecycle by subclasses. In the standard Xolotl clusters, this
	 * vector is filled in dissociateCluster that is called by
	 * createDissociationConnectivity.
	 */
	std::vector<ClusterPair> dissociatingPairs;
	std::vector<ClusterPair0> dissociatingPairs0;

	/**
	 * A vector of ClusterPairs that represent pairs of clusters that are emitted
	 * from the dissociation of this cluster. This vector should be populated early
	 * in the cluster's lifecycle by subclasses. In the standard Xolotl clusters,
	 * this vector is filled in emitClusters that is called by
	 * createDissociationConnectivity.
	 */
	std::vector<ClusterPair> emissionPairs;
	std::vector<ClusterPair0> emissionPairs0;

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	PSICluster() = delete;

	/**
	 * Construct a PSICluster.
	 *
	 * @param registry The performance handler registry
	 */
	PSICluster(PSIClusterReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
			const std::string& _name = "PSICluster");

	/**
	 * Copy constructor.
	 */
	PSICluster(PSICluster &other) :
			Reactant(other),
            indexList(other.indexList),
            psDim(other.psDim) {
	}

	/**
	 * The destructor
	 */
	virtual ~PSICluster() {
	}

	/**
	 * Update reactant using other reactants in its network.
	 */
	virtual void updateFromNetwork() override;

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	void resultFrom(ProductionReaction& reaction, int a[4] = defaultInit,
			int b[4] = defaultInit) override;

	/**
	 * Note that we result from the given reaction involving a super cluster.
	 * Assumes the reaction is already in the network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param prInfos Production reaction parameters.
	 */
	void resultFrom(ProductionReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param product The cluster created by the reaction.
	 *
	 */
	void resultFrom(ProductionReaction& reaction, IReactant& product) override;

	/**
	 * Note that we result from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param coef The cooresponding coefficient
	 */
	void resultFrom(ProductionReaction& reaction, double *coef) override;

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param a Number that can be used by daughter classes.
	 */
	void participateIn(ProductionReaction& reaction, int a[4] = defaultInit)
			override;

	/**
	 * Note that we combine with another cluster in a production reaction
	 * involving a super cluster.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param prInfos Production reaction parameters.
	 */
	void participateIn(ProductionReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we combine with another cluster in a production reaction
	 * involving a super cluster.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param product The cluster created by the reaction.
	 */
	void participateIn(ProductionReaction& reaction, IReactant& product)
			override;

	/**
	 * Note that we combine with another cluster in a production reaction.
	 * Assumes that the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster takes part.
	 * @param coef Number that can be used by daughter classes.
	 */
	void participateIn(ProductionReaction& reaction, double *coef) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param a Number that can be used by daughter classes.
	 * @param b Number that can be used by daughter classes.
	 */
	void participateIn(DissociationReaction& reaction, int a[4] = defaultInit,
			int b[4] = defaultInit) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction
	 * involving a super cluster.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param prInfos Production reaction parameters.
	 */
	void participateIn(DissociationReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction
	 * involving a super cluster.
	 * Assumes the reaction is already inour network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param disso The dissociating cluster.
	 */
	void participateIn(DissociationReaction& reaction, IReactant& disso)
			override;

	/**
	 * Note that we combine with another cluster in a dissociation reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction creating this cluster.
	 * @param coef Number that can be used by daughter classes.
	 */
	void participateIn(DissociationReaction& reaction, double *coef) override;

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param a Number that can be used by daughter classes.
	 */
	void emitFrom(DissociationReaction& reaction, int a[4] = defaultInit)
			override;

	/**
	 * Note that we emit from the given reaction involving a super cluster.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param prInfos Production reaction parameters.
	 */
	void emitFrom(DissociationReaction& reaction,
			const std::vector<PendingProductionReactionInfo>& prInfos) override;

	/**
	 * Note that we emit from the given reaction involving a super cluster.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param disso The dissociating cluster.
	 */
	void emitFrom(DissociationReaction& reaction, IReactant& disso) override;

	/**
	 * Note that we emit from the given reaction.
	 * Assumes the reaction is already in our network.
	 *
	 * @param reaction The reaction where this cluster emits.
	 * @param coef Number that can be used by daughter classes.
	 */
	void emitFrom(DissociationReaction& reaction, double *coef) override;

	/**
	 * This operation returns the connectivity array for this cluster for
	 * forward reactions. An entry with value one means that this cluster
	 * and the cluster with id = index + 1 are connected.
	 * 
	 * @return The connectivity array for "forward" (non-dissociating)
	 * reactions
	 */
	virtual std::vector<int> getReactionConnectivity() const;

	/**
	 * This operation returns the connectivity array for this cluster for
	 * forward reactions. An entry with value one means that this cluster
	 * and the cluster with id = index + 1 are connected.
	 * 
	 * @return The connectivity array for "forward" (non-dissociating)
	 * reactions
	 */
	virtual std::vector<int> getDissociationConnectivity() const;

	/**
	 * This operation returns the first moment.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param axis The axis we are intersted in
	 * @return The moment
	 */
	virtual double getMoment(const double* __restrict concs, int axis) const {
		return 0.0;
	}

	/**
	 * This operation returns the distance to the mean.
	 *
	 * @param atom The number of atoms
	 * @param axis The axis we are intersted in
	 * @return The distance to the mean number of atoms in the group
	 */
	virtual double getDistance(int atom, int axis) const {
		return 0.0;
	}

	/**
	 * This operation returns the factor used for the moments.
	 *
	 * @param atom The number of atoms
	 * @param axis The axis we are intersted in
	 * @return The factor
	 */
	virtual double getFactor(int atom, int axis) const {
		return 0.0;
	}

	/**
	 * This operation works as getPartialDerivatives above, but instead of
	 * returning a vector that it creates it fills a vector that is passed to
	 * it by the caller. This allows the caller to optimize the amount of
	 * memory allocations to just one if they are accessing the partial
	 * derivatives many times.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] partials The vector that should be filled with the partial
     * derivatives for this reactant where index zero corresponds to the 
     * first reactant in the list returned by the ReactionNetwork::getAll()
     * operation. The size of the vector should be equal to 
     * ReactionNetwork::size().
	 */
	virtual void getPartialDerivatives(const double* __restrict concs, int i,
            std::vector<double> & partials) const override;

    virtual void computePartials0(const double* __restrict concs, int xi,
            std::vector<double>& partials) const {

        // Get the partial derivatives for each reaction type
        computeAllProdPartials0(concs, xi, partials);
        computeAllCombPartials0(concs, xi, partials);
        computeAllDissPartials0(xi, partials);
        computeAllEmitPartials0(xi, partials);
    }


	/**
	 * This operation computes the partial derivatives due to production
	 * reactions.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param[out] partials The vector into which the partial derivatives 
     * should be inserted. This vector should have a length equal to the 
     * size of the network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getProductionPartialDerivatives(const double* __restrict concs, int i,
			std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to combination
	 * reactions.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param[out] partials The vector into which the partial derivatives 
     * should be inserted. This vector should have a length equal to the 
     * size of the network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getCombinationPartialDerivatives(const double* __restrict concs, int i,
			std::vector<double> & partials) const;

	/**
	 * This operation computes the partial derivatives due to dissociation of
	 * other clusters into this one.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getDissociationPartialDerivatives(
			std::vector<double> & partials, int i) const;

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	virtual void getEmissionPartialDerivatives(
			std::vector<double> & partials, int i) const;

	/**
	 * This operation reset the connectivity sets based on the information
	 * in the effective production and dissociation vectors.
	 */
	void resetConnectivities() override;

	/**
	 * This operation returns the sum of combination rate and emission rate
	 * (where this cluster is on the left side of the reaction) for this
	 * particular cluster.
	 * This is used to computed the desorption rate in the
	 * modified trap-mutation handler.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The position on the grid
	 * @return The rate
	 */
	double getLeftSideRate(const double* __restrict concs, int i) const override;

	/**
	 * This operation returns the vector of production reactions in which
	 * this cluster is involved, containing the id of the reactants, and
	 * the coefs[0][0]
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getProdVector() const override;

	/**
	 * This operation returns the vector of combination reactions in which
	 * this cluster is involved, containing the id of the other reactants, and
	 * the coefs[0]
	 *
	 * @return The vector of combinations
	 */
	virtual std::vector<std::vector<double> > getCombVector() const override;

	/**
	 * This operation returns the vector of dissociation reactions in which
	 * this cluster is involved, containing the id of the emitting reactants, and
	 * the coefs[0][0]
	 *
	 * @return The vector of dissociations
	 */
	virtual std::vector<std::vector<double> > getDissoVector() const override;

	/**
	 * This operation returns the vector of emission reactions in which
	 * this cluster is involved, containing
	 * the coefs[0][0]
	 *
	 * @return The vector of productions
	 */
	virtual std::vector<std::vector<double> > getEmitVector() const override;

	/**
	 * This operation returns a list that represents the connectivity
	 * between this cluster and other clusters in the network.
	 * "Connectivity" indicates whether two clusters interact, via any
	 * mechanism, in an abstract sense (as if they were nodes connected by
	 * an edge on a network graph).
	 *
	 * @return An array of ones and zeros that indicate whether or not this
	 * cluster interacts via any mechanism with another cluster. A "1" at
	 * the i-th entry in this array indicates that the cluster interacts
	 * with the i-th cluster in the ReactionNetwork and a "0" indicates
	 * that it does not.
	 */
	std::vector<int> getConnectivity() const override;

	/**
	 * Access bounds on number of given atoms represented by this cluster.
	 *
	 * @ param axis The direction
	 */
	// TODO do we want to make this generic by taking a type parameter?
	const IntegerRange<IReactant::SizeType>& getBounds(int axis) const {
		return bounds[axis];
	}

	/**
	 * Tell reactant to output a representation of its reaction coefficients
	 * to the given output stream.
	 *
	 * @param os Output stream on which to output coefficients.
	 */
	virtual void outputCoefficientsTo(std::ostream& os) const override;


	/**
     * Compute total flux(es) of this reactant using current concentrations
     * into their respective locations in the output concentrations.
	 *
     * @param concs Current concentrations for desired grid point.
	 * @param xi The location on the grid in the depth direction
     * @param updatedConcs Updated concentrations for desired grid point.
	 */
    void computeTotalFluxes(const double* __restrict concs, int xi,
                            double* __restrict updatedConcs) const override {
        // Compute the total fluxes for reactions we participate in.
        auto flux = (psDim == 1) ?
            getTotalFluxHelper0<Reactant::Flux>(concs, xi) :
            getTotalFluxHelper<Reactant::Flux>(concs, xi);

        // Update our concentration in the output concentration array.
        addToConcentration(updatedConcs, flux.flux);
    }

    virtual void useZerothMomentSpecializations();
};

} /* end namespace xolotlCore */
#endif
