#pragma once

#include "double_comparing.h"
#include "matrix_data.h"

namespace Matrices {

template <internal::matrix_elem T>
class matrix_t final : public internal::matrix_data<T> {
    using internal::matrix_data<T>::n_;
    using internal::matrix_data<T>::access;

public:
    explicit matrix_t(std::size_t n, const T &val = T(0))
        : internal::matrix_data<T>(n, val)
    {}

    matrix_t(const std::initializer_list<T> &inpt)
        : internal::matrix_data<T>(inpt)
    {}

    template <std::forward_iterator FwdIt>
    matrix_t(FwdIt start, FwdIt fin) : internal::matrix_data<T>(start, fin)
    {}

private:
    int make_non_zero_row(std::size_t row_id) noexcept
    {
        if (!DblCmp::are_eq(*access(row_id, row_id), T(0)))
            return 1;
        for (std::size_t i = row_id + 1; i < n_; ++i)
        {
            if (!DblCmp::are_eq(*access(i, row_id), T(0)))
            {
                this->swap_rows(row_id, i);
                return -1;
            }
        }
        return 0;
    }

    int normalize_rows(std::size_t row_id) noexcept
    {
        assert(!DblCmp::are_eq(*access(row_id, row_id), T(0)));
        int transform_coeff = 1;
        for (std::size_t i = row_id + 1; i < n_; ++i)
        {
            if (DblCmp::are_eq(*access(i, row_id), T(0)))
                continue;

            if (DblCmp::are_geq(std::abs(*access(i, row_id)),
                                std::abs(*access(row_id, row_id))))
            {
                transform_coeff *= -1;
                this->swap_rows(i, row_id);
            }

            T norm_coeff = *access(i, row_id) / *access(row_id, row_id);
            this->sub_row(i, row_id, norm_coeff);
        }
        return transform_coeff;
    }

    T calculate_det_echelon() const
    {
        matrix_t<T> tmp{*this};
        int det_coeff = tmp.make_echelon_form();
        return det_coeff * tmp.trace();
    }

public:
    int make_echelon_form()
    {
        int transform_coeff = 1;
        for (std::size_t i = 0; i < n_; ++i)
        {
            transform_coeff *= make_non_zero_row(i);
            if (transform_coeff == 0)
                return transform_coeff;
            transform_coeff *= normalize_rows(i);
        }
        return transform_coeff;
    }

    T calculate_det() const
    {
        if (n_ == 0)
            return T(1);
        return calculate_det_echelon();
    }
};

template <internal::integral_elem T>
class matrix_t<T> final : public internal::matrix_data<T> {
    using internal::matrix_data<T>::n_;
    using internal::matrix_data<T>::access;

public:
    explicit matrix_t(std::size_t n, const T &val = T(0))
        : internal::matrix_data<T>(n, val)
    {}

    matrix_t(const std::initializer_list<T> &inpt)
        : internal::matrix_data<T>(inpt)
    {}

    template <std::forward_iterator FwdIt>
    matrix_t(FwdIt start, FwdIt fin) : internal::matrix_data<T>(start, fin)
    {}

private:
    T calculate_det_echelon() const
    {
        // Sadly, only cast to double.
        matrix_t<double> tmp{this->cbegin(), this->cend()};
        int det_coeff = tmp.make_echelon_form();
        return det_coeff * tmp.trace();
    }

public:
    T calculate_det() const
    {
        if (n_ == 0)
            return T(1);
        return calculate_det_echelon();
    }
};

} // namespace Matrices