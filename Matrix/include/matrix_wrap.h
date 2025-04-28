#pragma once

#include "matrix.h"

#include <type_traits>

namespace Matrices {

namespace internal {

template <typename T> concept wrap_elem = std::is_convertible_v<T, double>;
template <typename T> concept integral_elem = std::is_integral_v<T>;

} // namespace internal

template <typename T> class matrix_wrap_t final : public matrix_t<T> {
public:
    explicit matrix_wrap_t(std::size_t n, const T &val = T(0))
        : matrix_t<T>(n, val)
    {}

    matrix_wrap_t(const std::initializer_list<T> &inpt) : matrix_t<T>(inpt) {}

    template <std::forward_iterator FwdIt>
    matrix_wrap_t(FwdIt start, FwdIt fin) : matrix_t<T>(start, fin)
    {}
};

template <internal::integral_elem T>
class matrix_wrap_t<T> final : public matrix_t<double> {
public:
    T calculate_det() const
    {
        return static_cast<T>(matrix_t<double>::calculate_det());
    }

    matrix_wrap_t &operator=(const std::initializer_list<T> &inpt)
    {
        matrix_t<double>::operator=(inpt);
        return *this;
    }

    matrix_wrap_t &operator*=(const T &val)
    {
        matrix_t<double>::operator*=(val);
        return *this;
    }

public:
    explicit matrix_wrap_t(std::size_t n, const T &val = T(0))
        : matrix_t<double>(n, static_cast<double>(val))
    {}

    matrix_wrap_t(const std::initializer_list<T> &inpt)
        : matrix_t<double>(inpt.begin(), inpt.end())
    {}

    template <std::forward_iterator FwdIt>
    matrix_wrap_t(FwdIt start, FwdIt fin) : matrix_t<double>(start, fin)
    {}
};

} // namespace Matrices