#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H

// Includes
#include <set>
#include "IReactionNetwork.h"
#include <Constants.h>

namespace xolotlPerf {
class IHandlerRegistry;
class IEventCounter;
}

// Override the hash operation for the composition maps used by the
// PSIClusterReactionNetwork to store reactants.
namespace std {
template<>
class hash<std::map<std::string, int>> {
public:
	long operator()(const std::map<std::string, int>& composition) const {
		int bigNumber = 1e9;
		return (composition.at(xolotlCore::heType) * 10
				+ composition.at(xolotlCore::vType) * 200
				+ composition.at(xolotlCore::iType) * 3000) * bigNumber;
	}
};
}

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants
 *  (combinations of normal reactants). It also manages a set of properties
 *  that describe both.
 */
class ReactionNetwork : public IReactionNetwork {

protected:

    /**
     * A functor useful for identifying a set of reactants by their
     * composition from a container, e.g., when removing a collection
     * of doomed reactants from a vector of reactants.
     */
    class ReactantMatcher
    {
    private:
        /**
         * The canonical composition string representations of the reactants
         * we are to find.
         */
        std::set<std::string> compStrings;

    public:
        /**
         * Build a ReactantMatcher.
         * @param reactants The collection of reactants we are to recognize.
         */
        ReactantMatcher(const std::vector<IReactant*>& reactants)
        {
            for(auto reactant : reactants)
            {
                compStrings.insert(reactant->getCompositionString());
            }
        }

        /**
         * Determine if the given reactant is in our set.
         * @param testReactant The reactant to check.
         * @return true iff the reactant's composition's canonical string
         * representation is in our set.
         */
        bool operator()(const IReactant* testReactant) const
        {
            auto iter = compStrings.find(testReactant->getCompositionString());
            return (iter != compStrings.end());
        }

        /**
         * Determine if the given reactant is in our set.
         * @param testReactant The reactant to check.
         * @return true iff the reactant's composition's canonical string
         * representation is in our set.
         */
        bool operator()(const std::shared_ptr<IReactant> testReactant) const
        {
            return this->operator()(testReactant.get());
        }
    };



	/**
	 * The properties of this network. The exact configuration of the map is
	 * specified by the class that loaded the network.
	 */
	std::shared_ptr<std::map<std::string, std::string>> properties;

	/**
	 * The performance handler registry that will be used with
	 * this class.
	 */
	std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

	/**
	 * Counter for the number of times the network concentration is updated.
	 */
	std::shared_ptr<xolotlPerf::IEventCounter> concUpdateCounter;

	/**
	 * The list of all of the reactants in the network. This list is filled and
	 * maintained by the getAll() operation.
	 */
	std::shared_ptr< std::vector<IReactant *> > allReactants;

	/**
	 * The current temperature at which the network's clusters exist.
	 */
	double temperature;

	/**
	 * The size of the network. It is also used to set the id of new reactants
	 * that are added to the network.
	 */
	int networkSize;

	/**
	 * The names of the reactants supported by this network.
	 */
	std::vector<std::string> names;

	/**
	 * The names of the compound reactants supported by this network.
	 */
	std::vector<std::string> compoundNames;

	/**
	 * The default constructor. It initializes the properties map and reactants vector.
	 */
	ReactionNetwork();

public:

	/**
	 * The constructor that takes the performance handler registry.
	 * It initializes the properties map and reactants vector.
	 *
	 * @param registry The performance handler registry
	 */
	ReactionNetwork(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * The copy constructor.
	 *
	 * @param other The ReactionNetwork to copy
	 */
	ReactionNetwork(const ReactionNetwork &other);

	/**
	 * The destructor.
	 */
	virtual ~ReactionNetwork() {
	}

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants.
	 *
	 * @param temp The new temperature
	 */
	virtual void setTemperature(double temp);

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 *
	 * @return The temperature
	 */
	virtual double getTemperature() const;

	/**
	 * This operation returns a reactant with the given type and size if it
	 * exists in the network or null if not.
	 *
	 * @param type The type of the reactant
	 * @param size The size of the reactant
	 * @return A pointer to the reactant
	 */
	virtual IReactant * get(const std::string& type, const int size) const;

	/**
	 * This operation returns a compound reactant with the given type and size if it
	 * exists in the network or null if not.
	 *
	 * @param type The type of the compound reactant
	 * @param sizes An array containing the sizes of each piece of the reactant
	 * @return A pointer to the compound reactant
	 */
	virtual IReactant * getCompound(const std::string& type,
			const std::vector<int>& sizes) const;

	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 *
	 * @return The list of all of the reactants in the network
	 */
	virtual const std::shared_ptr<std::vector<IReactant *>> & getAll() const;

	/**
	 * This operation returns all reactants in the network with the given type.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 *
	 * @param type The reactant or compound reactant type
	 * @return The list of all of the reactants in the network or null if the
	 * type is invalid
	 */
	virtual std::vector<IReactant *> getAll(const std::string& type) const;

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 *
	 * @param reactant The reactant that should be added to the network
	 */
	virtual void add(std::shared_ptr<IReactant> reactant) {return;}

	/**
	 * This operation adds a super reactant to the network.
	 * Used with a grouping method.
	 *
	 * @param reactant The reactant that should be added to the network
	 */
	virtual void addSuper(std::shared_ptr<IReactant> reactant) {return;}

	/**
	 * This operation removes the reactants from the network.
	 *
	 * @param reactant The reactants that should be removed.
	 */
	virtual void removeReactants(const std::vector<IReactant*>& reactant) {return;}

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	virtual void reinitializeNetwork() {return;}

	/**
	 * This operation returns the names of the reactants in the network.
	 *
	 * @return A vector with one entry for each of the distinct reactant types
	 * in the network
	 */
	const std::vector<std::string> & getNames() const;

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 *
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network
	 */
	const std::vector<std::string> & getCompoundNames() const;

	/**
	 * This operation returns a map of the properties of this reaction network.
	 *
	 * @return The map of properties that has been configured for this
	 * ReactionNetwork.
	 */
	const std::map<std::string, std::string> & getProperties();

	/**
	 * This operation sets a property with the given key to the specified value
	 * for the network. ReactionNetworks may reserve the right to ignore this
	 * operation for special key types, most especially those that they manage
	 * on their own.
	 *
	 * @param key The key for the property
	 * @param value The value to which the key should be set
	 */
	virtual void setProperty(const std::string& key,
			const std::string& value) {return;}

	/**
	 * This operation returns the size or number of reactants in the network.
	 *
	 * @return The number of reactants in the network
	 */
	int size();

	/**
	 * This operation fills an array of doubles with the concentrations of all
	 * of the reactants in the network.
	 *
	 * @param concentrations The array that will be filled with the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. If the array is too small to hold the concentrations, SIGSEGV will
	 * be thrown.
	 */
	void fillConcentrationsArray(double * concentrations);

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	virtual void updateConcentrationsFromArray(double * concentrations);

	/**
	 * Request that all reactants in the network release their
	 * pointers to the network, to break cycles and allow the
	 * network and the clusters/reactant objects it contains to
	 * be destroyed gracefully.
	 *
	 * Should only be done when the network is no longer needed.
	 * Ideally, we would do this from the network's destructor.
	 * However, unless the reactants in the network release their
	 * shared_ptrs to the network, the reference count on the
	 * network's shared_ptr will never reach zero and the
	 * object it owns will never be destroyed.  (Hence, the
	 * reactant objects will also never be destroyed.)
	 */
	void askReactantsToReleaseNetwork();

};

}

#endif
