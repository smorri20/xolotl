#ifndef XCORE_REACTANTS_PSI_COEFFS_FINDER_H
#define XCORE_REACTANTS_PSI_COEFFS_FINDER_H

namespace xolotlCore {

template<typename T>
struct CoeffsFinder {

    using ValueType = T;

    // Info about the reactions being defined.
    const std::vector<PendingProductionReactionInfo>& prInfos;

    // Construct a coefficients finder.
    CoeffsFinder() = delete;
    CoeffsFinder(const CoeffsFinder& other) = default;
    CoeffsFinder(const std::vector<PendingProductionReactionInfo>& _prInfos)
      : prInfos(_prInfos)
    { }

    // Initialize our running value for a reduction.
    KOKKOS_INLINE_FUNCTION
    void init(ValueType& running) const {
        running.fill(0);
    }

    // Combine two running values.
    KOKKOS_INLINE_FUNCTION
    void join(volatile ValueType& dest,
                const volatile ValueType& src) const {

        dest += src;
    }
};

} // namespace xolotlCore

#endif // XCORE_REACTANTS_PSI_COEFFS_FINDER_H
