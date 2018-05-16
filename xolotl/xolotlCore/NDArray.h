#ifndef XCORE_NDARRAY_H
#define XCORE_NDARRAY_H

#include <iostream>
#include <sstream>
#include <iterator>
#include "Kokkos_Core.hpp"

namespace xolotlCore {

// Forward declaration, apparently needed for the rest to
// compile successfully with some compilers (e.g., clang++ 6).
template<typename T, uint32_t N0, uint32_t... Ns>
class Array;


// An N-dimensional array of Ts.
// Specialization for 2+ dimensions.
template<typename T, uint32_t N0, uint32_t N1, uint32_t... Ns>
class Array<T, N0, N1, Ns...>
{
public:
    using value_type = Array<T, N1, Ns...>;
    using size_type = std::size_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using iterator = pointer;
    using const_iterator = const_pointer;

private:
    value_type data[N0];

public:    
    Array(void) = default;
    Array(const Array<T, N0, N1, Ns...>& other) = default;

    KOKKOS_INLINE_FUNCTION
    void fill(const T& val)
    {
        for(auto i = 0; i < N0; ++i)
        {
            data[i].fill(val);
        }
    }

    iterator begin(void) noexcept
    {
        return data;
    }

    const_iterator begin(void) const noexcept
    {
        return data;
    }

    iterator end(void) noexcept
    {
        return data + N0;
    }

    const_iterator end(void) const noexcept
    {
        return data + N0;
    }

    constexpr size_type max_size(void) const noexcept
    {
        return N0;
    }

    constexpr size_type size(void) const noexcept
    {
        return N0;
    }

    constexpr bool empty(void) const noexcept
    {
        return (N0 == 0);
    }

    KOKKOS_INLINE_FUNCTION
    volatile Array<T, N1, Ns...>& operator[](size_t i) volatile
    {
        return data[i];
    }

    KOKKOS_INLINE_FUNCTION
    const volatile Array<T, N1, Ns...>& operator[](size_t i) const volatile
    {
        return data[i];
    }

    KOKKOS_INLINE_FUNCTION
    volatile Array<T, N0, N1, Ns...>& operator+=(const volatile Array<T, N0, N1, Ns...>& other) volatile
    {
        for(auto i = 0; i < N0; ++i)
        {
            data[i] += other[i];
        }
        return *this;
    }

    void DumpTo(std::ostream& os, size_t level = 0)
    {
        std::ostringstream indentStr;
        for(auto j = 0; j < level; ++j)
        {
            indentStr << ' ';
        }

        os << indentStr.str() << "[\n";
        for(auto i = 0; i < N0; ++i)
        {
            data[i].DumpTo(os, level + 1);
        }
        os << indentStr.str() << "]\n";
    }
};

template<typename T, uint32_t N0, uint32_t N1, uint32_t...Ns>
Array<T, N0, N1, Ns...>
operator+(const Array<T, N0, N1, Ns...>& a, const Array<T, N0, N1, Ns...>& b)
{
    Array<T, N0, N1, Ns...> ret = a;
    ret += b;
    return ret;
}


// An N-dimensional array of type T.
// Specialization for one dimension.
template<typename T, uint32_t N0>
class Array<T, N0>
{
public:
    using value_type = T;
    using size_type = std::size_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using iterator = pointer;
    using const_iterator = const_pointer;

private:
    value_type data[N0];

public:    
    Array(void) = default;
    Array(const Array<T, N0>& other) = default;

    KOKKOS_INLINE_FUNCTION
    void fill(const T& val) volatile
    {
        for(auto i = 0; i < N0; ++i)
        {
            data[i] = val;
        }
    }

    iterator begin(void) noexcept
    {
        return data;
    }

    const_iterator begin(void) const noexcept
    {
        return data;
    }

    iterator end(void) noexcept
    {
        return data + N0;
    }

    const_iterator end(void) const noexcept
    {
        return data + N0;
    }

    constexpr size_type max_size(void) const noexcept
    {
        return N0;
    }

    constexpr size_type size(void) const noexcept
    {
        return N0;
    }

    constexpr bool empty(void) const noexcept
    {
        return (N0 == 0);
    }

    KOKKOS_INLINE_FUNCTION
    volatile value_type& operator[](size_t i) volatile
    {
        return data[i];
    }

    KOKKOS_INLINE_FUNCTION
    const volatile value_type& operator[](size_t i) const volatile
    {
        return data[i];
    }

    KOKKOS_INLINE_FUNCTION
    volatile Array<T, N0>& operator+=(const volatile Array<T, N0>& other) volatile
    {
        for(auto i = 0; i < N0; ++i)
        {
            data[i] += other[i];
        }
        return *this;
    }

    void DumpTo(std::ostream& os, size_t level = 0)
    {
        for(auto j = 0; j < level; ++j)
        {
            os << ' ';
        }
        std::copy(this->begin(), this->end(),
                std::ostream_iterator<T>(os, " "));
        os << std::endl;
    }
};

template<typename T, uint32_t N0>
Array<T, N0>
operator+(const Array<T, N0>& a, const Array<T, N0>& b)
{
    Array<T, N0> ret = a;
    ret += b;
    return ret;
}

} // namespace xolotlCore

#endif // XCORE_NDARRAY_H
