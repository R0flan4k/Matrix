#pragma once

#include "double_comparing.h"
#include "matrix_exceptions.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iostream>
#include <iterator>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

using DblCmp::are_eq;
using DblCmp::are_geq;
using DblCmp::is_zero;

namespace Matrices {

struct matrix_dumper {
    template <typename MatrT>
    static void dump(const MatrT &m, const std::string &descr = "")
    {
        std::cout << "\t" << descr << std::endl;
        auto it = m.begin();
        for (std::size_t i = 0; i < m.rank(); ++i)
        {
            for (std::size_t j = 0; j < m.rank(); ++j, ++it)
                std::cout << *it << "\t";
            std::cout << std::endl;
        }
    }
};

template <typename T> concept matrix_elem = requires(T t)
{
    requires std::is_floating_point_v<T>;
    requires std::is_nothrow_constructible_v<T>;
    requires std::is_nothrow_move_constructible_v<T>;
};

template <matrix_elem T> class matrix_buff {
protected:
    std::size_t sz_;
    std::size_t used_;
    T *data_;

protected:
    static void construct(T *p, const T &val) { new (p) T(val); }
    static void construct(T *p, T &&val) { new (p) T(std::move(val)); }
    static void destroy(T *p) { p->~T(); }

    explicit matrix_buff(std::size_t sz)
        : sz_(sz), used_(0),
          data_((sz_ == 0) ? nullptr
                           : static_cast<T *>(::operator new(sizeof(T) * sz_)))
    {}

    explicit matrix_buff(const matrix_buff &m) : matrix_buff(m.sz_)
    {
        assert(m.used_ == m.sz_);
        for (std::size_t i = 0; i < sz_; ++i)
        {
            construct(data_ + i, m.data_[i]);
            ++used_;
        }
    }

    matrix_buff(matrix_buff &&m) noexcept
        : sz_(m.sz_), used_(m.used_), data_(m.data_)
    {}

    matrix_buff &operator=(const matrix_buff &m)
    {
        if (&m == this)
            return *this;
        matrix_buff<T> tmp(m);
        std::swap(data_, tmp.data_);
        std::swap(sz_, tmp.sz_);
        std::swap(used_, tmp.used_);
        return *this;
    }

    matrix_buff &operator=(matrix_buff &&m) noexcept
    {
        if (&m == this)
            return *this;
        std::swap(data_, m.data_);
        std::swap(sz_, m.sz_);
        std::swap(used_, m.used_);
        return *this;
    }

    virtual ~matrix_buff()
    {
        for (std::size_t i = 0; i < used_; ++i)
            destroy(data_ + i);
        ::operator delete(data_, sizeof(T) * sz_);
    }
};

template <matrix_elem T> class matrix_t : private matrix_buff<T> {
protected:
    using matrix_buff<T>::data_;
    using matrix_buff<T>::sz_;
    using matrix_buff<T>::used_;
    size_t n_;

public:
    typedef T *iterator;
    typedef const T *const_iterator;

protected:
    T *access(std::size_t row, std::size_t col) const
    {
        assert(row < n_ && col < n_);
        return data_ + n_ * row + col;
    }

public:
    std::size_t rank() const { return n_; }

    explicit matrix_t(std::size_t n, const T &val = T(0))
        : matrix_buff<T>(n * n), n_(n)
    {
        for (std::size_t i = 0; i < sz_; ++i)
        {
            this->construct(data_ + i, val);
            ++used_;
        }
        assert(used_ == sz_);
    }

    matrix_t(const std::initializer_list<T> &inpt)
        : matrix_buff<T>(inpt.size()), n_(std::sqrt(inpt.size()))
    {
        std::size_t i = 0;
        for (auto start = inpt.begin(), fin = inpt.end(); start != fin;
             ++start, ++i)
        {
            this->construct(data_ + i, *start);
            ++used_;
        }
        assert(used_ == sz_);
    }

    template <std::random_access_iterator RandIt>
    explicit matrix_t(RandIt start, RandIt fin)
        : matrix_buff<T>(std::distance(start, fin)),
          n_(std::sqrt(std::distance(start, fin)))
    {
        for (std::size_t i = 0, dist = std::distance(start, fin); i < dist;
             ++i, ++start)
        {
            this->construct(data_ + i, *start);
            ++used_;
        }
        assert(used_ == sz_);
    }

    virtual ~matrix_t() = default;

    matrix_t &operator=(const std::initializer_list<T> &inpt)
    {
        if (inpt.size() != sz_)
            throw MatrExcepts::wrong_init_list(
                "Wrong size of initializer list.");
        std::move(inpt.cbegin(), inpt.cend(), data_) assert(used_ == sz_);
        return *this;
    }

    matrix_t &operator*(const T &val) noexcept
    {
        for (std::size_t i = 0; i < sz_; ++i)
            data_[i] *= val;
        return *this;
    }

    void swap_rows(std::size_t lhs, std::size_t rhs) noexcept
    {
        if (lhs == rhs)
            return;
        for (std::size_t i = 0; i < n_; ++i)
            std::swap(*access(lhs, i), *access(rhs, i));
    }

    void mul_row(std::size_t row_id, const T &val) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(row_id, i) *= val;
    }

    void div_row(std::size_t row_id, const T &val) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(row_id, i) /= val;
    }

    void add_row(std::size_t lhs, std::size_t rhs, const T &val = T(1)) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(lhs, i) += *access(rhs, i) * val;
    }

    void sub_row(std::size_t lhs, std::size_t rhs, const T &val = T(1)) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(lhs, i) -= *access(rhs, i) * val;
    }

    template <std::random_access_iterator RandomIt>
    void set_col(std::size_t col_id, RandomIt start, RandomIt fin) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i, ++start)
            *access(i, col_id) = *start;
    }

protected:
    int make_non_zero_row(std::size_t row_id) noexcept
    {
        if (!are_eq(*access(row_id, row_id), T(0)))
            return 1;
        for (std::size_t i = row_id + 1; i < n_; ++i)
        {
            if (!are_eq(*access(i, row_id), T(0)))
            {
                swap_rows(row_id, i);
                return -1;
            }
        }
        return 0;
    }

    int normalize_rows(std::size_t row_id) noexcept
    {
        assert(!are_eq(*access(row_id, row_id), T(0)));
        int transform_coeff = 1;
        for (std::size_t i = row_id + 1; i < n_; ++i)
        {
            if (are_eq(*access(i, row_id), T(0)))
                continue;

            if (are_geq(std::abs(*access(i, row_id)),
                        std::abs(*access(row_id, row_id))))
            {
                transform_coeff *= -1;
                swap_rows(i, row_id);
            }

            T norm_coeff = *access(i, row_id) / *access(row_id, row_id);
            sub_row(i, row_id, norm_coeff);
        }
        return transform_coeff;
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

    T diag_multiplication() const noexcept
    {
        T res = T(1);
        for (std::size_t i = 0; i < n_; ++i)
            res *= *access(i, i);
        return res;
    }

    iterator begin() noexcept { return data_; }
    const_iterator begin() const noexcept { return data_; }
    iterator end() noexcept { return data_ + sz_; }
    const_iterator end() const noexcept { return data_ + sz_; }

protected:
    class row_t final {
        T *base_;

    public:
        row_t(T *base) noexcept : base_(base) {}

        const T *base() const { return base_; }

        T operator[](std::size_t off) const { return *(base_ + off); }

        T &operator[](std::size_t off) { return *(base_ + off); }
    };

public:
    row_t operator[](std::size_t off) const { return row_t{access(off, 0)}; }
};

template <matrix_elem T> class const_matrix_t : public matrix_t<T> {
    using matrix_t<T>::data_;
    using matrix_t<T>::sz_;
    using matrix_t<T>::used_;
    using matrix_t<T>::n_;
    std::optional<T> det_;

public:
    typedef matrix_t<T>::iterator iterator;
    typedef matrix_t<T>::const_iterator const_iterator;

public:
    explicit const_matrix_t(std::size_t n, const T &val = T(0))
        : matrix_t<T>(n, val)
    {}

    const_matrix_t(const std::initializer_list<T> &inpt) : matrix_t<T>(inpt) {}

    template <std::random_access_iterator RandIt>
    explicit const_matrix_t(RandIt start, RandIt fin) : matrix_t<T>(start, fin)
    {}

    const_matrix_t(const matrix_t<T> &matr) : matrix_t<T>(matr) {}

    const_matrix_t(matrix_t<T> &&matr) noexcept : matrix_t<T>(matr) {}

    virtual ~const_matrix_t() = default;

    const T &det() const
    {
        if (!det_.has_value())
            throw MatrExcepts::no_det("The determinant isn't calculated.");
        return det_.value();
    }

    matrix_t<T> &operator=(const std::initializer_list<T> &inpt) = delete;
    matrix_t<T> &operator*() = delete;
    void swap_rows() = delete;
    void mul_row() = delete;
    void add_row() = delete;
    void sub_row() = delete;
    void set_col() = delete;
    iterator begin() = delete;
    iterator end() = delete;

private:
    void calculate_det_echelon()
    {
        const_matrix_t<T> tmp{*this};
        int det_coeff = tmp.make_echelon_form();
        det_ = det_coeff * tmp.diag_multiplication();
    }

public:
    const T &calculate_det()
    {
        if (det_.has_value())
            return det_.value();
        if (n_ == 0)
            throw MatrExcepts::no_det(
                "Can't calculate empty matrix determinant.");
        calculate_det_echelon();
        return det_.value();
    }
};
} // namespace Matrices