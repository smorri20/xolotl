#ifndef XCORE_REACTANTS_PSI_SUPER_COEFFS_FINDER_H
#define XCORE_REACTANTS_PSI_SUPER_COEFFS_FINDER_H

#include "Kokkos_Core.hpp"
#include "CoeffsFinder.h"

namespace xolotlCore {

template<typename T>
struct SuperCoeffsFinder : public CoeffsFinder<T> {

    double numHe;
    double numV;
    double dispersionHe;
    double dispersionV;

    // Construct a coefficients finder.
    SuperCoeffsFinder() = delete;
    SuperCoeffsFinder(const SuperCoeffsFinder& other) = default;
    SuperCoeffsFinder(double _numHe, double _numV,
                    double _dispersionHe, double _dispersionV,
                    const std::vector<PendingProductionReactionInfo>& _prInfos)
      : CoeffsFinder<T>(_prInfos),
        numHe(_numHe),
        numV(_numV),
        dispersionHe(_dispersionHe),
        dispersionV(_dispersionV)
    { }
};


struct SuperProductCoeffsFinder : 
            public SuperCoeffsFinder<Array<double, 3, 3, 3>> {

    // Reactants involved in the reaction.
    PSICluster& r1;
    PSICluster& r2;

    // Construct a coefficients finder.
    SuperProductCoeffsFinder() = delete;
    SuperProductCoeffsFinder(const SuperProductCoeffsFinder& other) = default;
    SuperProductCoeffsFinder(PSICluster& _r1, PSICluster& _r2,
                    double _numHe, double _numV,
                    double _dispersionHe, double _dispersionV,
                    const std::vector<PendingProductionReactionInfo>& _prInfos)
      : SuperCoeffsFinder<Array<double, 3, 3, 3>>(_numHe, _numV,
                                                _dispersionHe, _dispersionV,
                                                _prInfos),
        r1(_r1),
        r2(_r2)
    { }


    KOKKOS_INLINE_FUNCTION
    void operator()(size_t idx, ValueType& running) const {

        // Use names corresponding to those in single version.
        auto const& currPRI = prInfos[idx];
        int a = currPRI.numHe;
        int b = currPRI.numV;
        int c = currPRI.i;
        int d = currPRI.j;

        double firstHeDistance = 0.0;
        double firstVDistance = 0.0;
        double secondHeDistance = 0.0;
        double secondVDistance = 0.0;
        if (r1.getType() == ReactantType::PSISuper) {
            firstHeDistance = r1.getHeDistance(c);
            firstVDistance = r1.getVDistance(d);
        }
        if (r2.getType() == ReactantType::PSISuper) {
            secondHeDistance = r2.getHeDistance(c);
            secondVDistance = r2.getVDistance(d);
        }
        double heFactor = (double) (a - numHe) / dispersionHe;
        double vFactor = (double) (b - numV) / dispersionV;
        // First is A, second is B, in A + B -> this
        running[0][0][0] += 1.0;
        running[0][0][1] += heFactor;
        running[0][0][2] += vFactor;
        running[1][0][0] += firstHeDistance;
        running[1][0][1] += firstHeDistance * heFactor;
        running[1][0][2] += firstHeDistance * vFactor;
        running[2][0][0] += firstVDistance;
        running[2][0][1] += firstVDistance * heFactor;
        running[2][0][2] += firstVDistance * vFactor;
        running[0][1][0] += secondHeDistance;
        running[0][1][1] += secondHeDistance * heFactor;
        running[0][1][2] += secondHeDistance * vFactor;
        running[0][2][0] += secondVDistance;
        running[0][2][1] += secondVDistance * heFactor;
        running[0][2][2] += secondVDistance * vFactor;
        running[1][1][0] += firstHeDistance * secondHeDistance;
        running[1][1][1] += firstHeDistance * secondHeDistance * heFactor;
        running[1][1][2] += firstHeDistance * secondHeDistance * vFactor;
        running[1][2][0] += firstHeDistance * secondVDistance;
        running[1][2][1] += firstHeDistance * secondVDistance * heFactor;
        running[1][2][2] += firstHeDistance * secondVDistance * vFactor;
        running[2][1][0] += firstVDistance * secondHeDistance;
        running[2][1][1] += firstVDistance * secondHeDistance * heFactor;
        running[2][1][2] += firstVDistance * secondHeDistance * vFactor;
        running[2][2][0] += firstVDistance * secondVDistance;
        running[2][2][1] += firstVDistance * secondVDistance * heFactor;
        running[2][2][2] += firstVDistance * secondVDistance * vFactor;
    }
};



struct SuperCombiningCoeffsFinder : 
            public SuperCoeffsFinder<Array<double, 3, 3, 3>> {

    // Reactants involved in the reaction.
    PSICluster& reactant;

    // Construct a coefficients finder.
    SuperCombiningCoeffsFinder() = delete;
    SuperCombiningCoeffsFinder(const SuperCombiningCoeffsFinder& other) = default;
    SuperCombiningCoeffsFinder(PSICluster& _reactant,
                    double _numHe, double _numV,
                    double _dispersionHe, double _dispersionV,
                    const std::vector<PendingProductionReactionInfo>& _prInfos)
      : SuperCoeffsFinder<Array<double, 3, 3, 3>>(_numHe, _numV,
                                                _dispersionHe, _dispersionV,
                                                _prInfos),
        reactant(_reactant)
    { }


    KOKKOS_INLINE_FUNCTION
    void operator()(size_t idx, ValueType& running) const {

        // Use names corresponding to those in single version.
        auto const& currPRI = prInfos[idx];
        int a = currPRI.i;
        int b = currPRI.j;

        double heDistance = reactant.getHeDistance(a);
        double heFactor = (double) (a - numHe) / dispersionHe;
        double vDistance = reactant.getVDistance(b);
        double vFactor = (double) (b - numV) / dispersionV;
        // This is A, itBis is B, in A + B -> C
        running[0][0][0] += 1.0;
        running[0][0][1] += heFactor;
        running[0][0][2] += vFactor;
        running[1][0][0] += heDistance;
        running[1][0][1] += heDistance * heFactor;
        running[1][0][2] += heDistance * vFactor;
        running[2][0][0] += vDistance;
        running[2][0][1] += vDistance * heFactor;
        running[2][0][2] += vDistance * vFactor;
    }
};



struct SuperDissociatingCoeffsFinder : 
            public SuperCoeffsFinder<Array<double, 3, 3>> {

    // Reactants involved in the reaction.
    IReactant& reactant;

    // Construct a coefficients finder.
    SuperDissociatingCoeffsFinder() = delete;
    SuperDissociatingCoeffsFinder(const SuperDissociatingCoeffsFinder& other) = default;
    SuperDissociatingCoeffsFinder(IReactant& _reactant,
                    double _numHe, double _numV,
                    double _dispersionHe, double _dispersionV,
                    const std::vector<PendingProductionReactionInfo>& _prInfos)
      : SuperCoeffsFinder<Array<double, 3, 3>>(_numHe, _numV,
                                                _dispersionHe, _dispersionV,
                                                _prInfos),
        reactant(_reactant)
    { }


    KOKKOS_INLINE_FUNCTION
    void operator()(size_t idx, ValueType& running) const {

        // Use names corresponding to the single-item version.
        auto const& currPRI = prInfos[idx];
        int a = currPRI.numHe;
        int b = currPRI.numV;
        int c = currPRI.i;
        int d = currPRI.j;

        double firstHeDistance = 0.0;
        double firstVDistance = 0.0;

        if (reactant.getType() == ReactantType::PSISuper) {
            auto const& super = static_cast<PSICluster const&>(reactant);
            firstHeDistance = super.getHeDistance(a);
            firstVDistance = super.getVDistance(b);
        }
        double heFactor = (double) (c - numHe) / dispersionHe;
        double vFactor = (double) (d - numV) / dispersionV;

        // A is the dissociating cluster
        running[0][0] += 1.0;
        running[0][1] += heFactor;
        running[0][2] += vFactor;
        running[1][0] += firstHeDistance;
        running[1][1] += firstHeDistance * heFactor;
        running[1][2] += firstHeDistance * vFactor;
        running[2][0] += firstVDistance;
        running[2][1] += firstVDistance * heFactor;
        running[2][2] += firstVDistance * vFactor;
    }
};



struct SuperEmittingCoeffsFinder : 
            public SuperCoeffsFinder<Array<double, 3, 3>> {

    // Reactants involved in the reaction.
    PSICluster& reactant;

    // Construct a coefficients finder.
    SuperEmittingCoeffsFinder() = delete;
    SuperEmittingCoeffsFinder(const SuperEmittingCoeffsFinder& other) = default;
    SuperEmittingCoeffsFinder(PSICluster& _reactant,
                    double _numHe, double _numV,
                    double _dispersionHe, double _dispersionV,
                    const std::vector<PendingProductionReactionInfo>& _prInfos)
      : SuperCoeffsFinder<Array<double, 3, 3>>(_numHe, _numV,
                                                _dispersionHe, _dispersionV,
                                                _prInfos),
        reactant(_reactant)
    { }


    KOKKOS_INLINE_FUNCTION
    void operator()(size_t idx, ValueType& running) const {

        // Use same names as used in single version.
        auto const& currPRI = prInfos[idx];
        int a = currPRI.numHe;
        int b = currPRI.numV;

        double heDistance = reactant.getHeDistance(a);
        double heFactor = (double) (a - numHe) / dispersionHe;
        double vDistance = reactant.getVDistance(b);
        double vFactor = (double) (b - numV) / dispersionV;
        // A is the dissociating cluster
        running[0][0] += 1.0;
        running[0][1] += heFactor;
        running[0][2] += vFactor;
        running[1][0] += heDistance;
        running[1][1] += heDistance * heFactor;
        running[1][2] += heDistance * vFactor;
        running[2][0] += vDistance;
        running[2][1] += vDistance * heFactor;
        running[2][2] += vDistance * vFactor;
    }
};

} // namespace xolotlCore

#endif // XCORE_REACTANTS_PSI_SUPER_COEFFS_FINDER_H
