#pragma once

#include "double_comparing.h"
#include "matrix_exceptions.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

namespace Matrices {

namespace internal {

template <typename T> concept matrix_like = requires(T t)
{
    {
        t.cbegin()
    }
    ->std::convertible_to<typename T::const_iterator>;
    {
        t.begin()
    }
    ->std::convertible_to<typename T::iterator>;
    {
        t.cend()
    }
    ->std::convertible_to<typename T::const_iterator>;
    {
        t.end()
    }
    ->std::convertible_to<typename T::iterator>;
    {
        t.rank()
    }
    ->std::same_as<std::size_t>;
};

template <typename T> concept matrix_elem = requires(T t)
{
    requires std::is_floating_point_v<T>;
    requires std::is_nothrow_constructible_v<T>;
    requires std::is_nothrow_move_constructible_v<T>;
};

template <class U, class T> concept t_like = std::is_convertible<U, T>::value;

template <typename T> concept pointer = std::is_pointer_v<T>;

template <internal::matrix_elem T> class matrix_buff {
protected:
    std::size_t sz_;
    std::size_t used_;
    T *data_;

protected:
    template <internal::t_like<T> U, typename K>
    static void construct(K *p, U &&val)
    {
        new (p) T(std::forward<U>(val));
    }

    static void destroy(T *p) { p->~T(); }

    void
    swap(matrix_buff &m) noexcept(std::is_nothrow_move_constructible_v<T>
                                      &&std::is_nothrow_move_assignable_v<T>)
    {
        std::swap(data_, m.data_);
        std::swap(sz_, m.sz_);
        std::swap(used_, m.used_);
    }

    explicit matrix_buff(std::size_t sz)
        : sz_(sz), used_(0),
          data_((sz_ == 0) ? nullptr
                           : static_cast<T *>(::operator new(sizeof(T) * sz_)))
    {}

    matrix_buff(const matrix_buff &m) : matrix_buff(m.sz_)
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
        swap(m);
        return *this;
    }

    matrix_buff &operator=(matrix_buff &&m) noexcept
    {
        if (&m == this)
            return *this;
        swap(m);
        return *this;
    }

    virtual ~matrix_buff()
    {
        for (std::size_t i = 0; i < used_; ++i)
            destroy(data_ + i);
        ::operator delete(data_, sizeof(T) * sz_);
    }
};

template <typename T, typename Distance = int, typename Pointer = T *,
          typename Reference = T &>
class random_access_iterator final {
    Pointer ptr_;

public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = Distance;
    using value_type = T;
    using reference = Reference;
    using pointer = Pointer;
    using size_type = std::size_t;

    random_access_iterator() noexcept : ptr_(Pointer()) {}
    random_access_iterator(const Pointer &ptr) noexcept : ptr_(ptr) {}

    random_access_iterator<T, Distance, Pointer, Reference> &
    operator+=(size_type n)
    {
        ptr_ += n;
        return *this;
    }
    random_access_iterator<T, Distance, Pointer, Reference> &
    operator-=(size_type n)
    {
        ptr_ -= n;
        return *this;
    }
    random_access_iterator<T, Distance, Pointer, Reference> operator++(int)
    {
        return random_access_iterator<T, Distance, Pointer, Reference>{ptr_++};
    }
    random_access_iterator<T, Distance, Pointer, Reference> &operator++()
    {
        return *this += 1;
    }
    random_access_iterator<T, Distance, Pointer, Reference> operator--(int)
    {
        return random_access_iterator<T, Distance, Pointer, Reference>{ptr_--};
    }
    random_access_iterator<T, Distance, Pointer, Reference> &operator--()
    {
        return *this -= 1;
    }
    reference operator*() const { return *ptr_; }
    pointer operator->() const { return ptr_; }

    friend difference_type operator-(
        const random_access_iterator<T, Distance, Pointer, Reference> &it1,
        const random_access_iterator<T, Distance, Pointer, Reference> &it2)
    {
        return it1.ptr_ - it2.ptr_;
    }

    friend bool operator==(
        const random_access_iterator<T, Distance, Pointer, Reference> &it1,
        const random_access_iterator<T, Distance, Pointer, Reference> &it2)
    {
        return it1.ptr_ == it2.ptr_;
    }

    reference operator[](size_type n) const { return ptr_[n]; }
};

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
bool operator!=(
    const random_access_iterator<T, Distance, Pointer, Reference> &it1,
    const random_access_iterator<T, Distance, Pointer, Reference> &it2)
{
    return !(it1 == it2);
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
bool operator>(
    const random_access_iterator<T, Distance, Pointer, Reference> &it1,
    const random_access_iterator<T, Distance, Pointer, Reference> &it2)
{
    return it1 - it2 > 0;
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
bool operator<(
    const random_access_iterator<T, Distance, Pointer, Reference> &it1,
    const random_access_iterator<T, Distance, Pointer, Reference> &it2)
{
    return it1 - it2 < 0;
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
bool operator>=(
    const random_access_iterator<T, Distance, Pointer, Reference> &it1,
    const random_access_iterator<T, Distance, Pointer, Reference> &it2)
{
    return it1 - it2 >= 0;
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
bool operator<=(
    const random_access_iterator<T, Distance, Pointer, Reference> &it1,
    const random_access_iterator<T, Distance, Pointer, Reference> &it2)
{
    return it1 - it2 <= 0;
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
random_access_iterator<T, Distance, Pointer, Reference>
operator+(random_access_iterator<T, Distance, Pointer, Reference> it,
          std::size_t n)
{
    return it += n;
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
random_access_iterator<T, Distance, Pointer, Reference>
operator+(std::size_t n,
          random_access_iterator<T, Distance, Pointer, Reference> it)
{
    return it += n;
}

template <typename T, typename Distance = std::size_t, typename Pointer = T *,
          typename Reference = T &>
random_access_iterator<T, Distance, Pointer, Reference>
operator-(random_access_iterator<T, Distance, Pointer, Reference> it,
          std::size_t n)
{
    return it -= n;
}

} // namespace internal

class matrix_dumper final {
    std::ostream &debug_stream_;

public:
    matrix_dumper(std::ostream &ds = std::cout) : debug_stream_(ds) {}

    template <internal::matrix_like MatrT>
    void dump(const MatrT &m, std::string_view descr = "") const
    {
        debug_stream_ << "\t" << descr << std::endl;
        auto it = m.cbegin();
        for (std::size_t i = 0; i < m.rank(); ++i)
        {
            for (std::size_t j = 0; j < m.rank(); ++j, ++it)
                debug_stream_ << *it << "\t";
            debug_stream_ << std::endl;
        }
    }
};

template <internal::matrix_elem T>
class matrix_t : private internal::matrix_buff<T> {
protected:
    using internal::matrix_buff<T>::data_;
    using internal::matrix_buff<T>::sz_;
    using internal::matrix_buff<T>::used_;
    size_t n_;

public:
    typedef internal::random_access_iterator<T> iterator;
    typedef internal::random_access_iterator<const T> const_iterator;

protected:
    T *access(std::size_t row, std::size_t col) const
    {
        assert(row < n_ + 1 && col < n_ + 1);
        return data_ + n_ * row + col;
    }

public:
    std::size_t rank() const { return n_; }

    explicit matrix_t(std::size_t n, const T &val = T(0))
        : internal::matrix_buff<T>(n * n), n_(n)
    {
        for (std::size_t i = 0; i < sz_; ++i)
        {
            this->construct(data_ + i, val);
            ++used_;
        }
        assert(used_ == sz_);
    }

    matrix_t(const std::initializer_list<T> &inpt)
        : internal::matrix_buff<T>(inpt.size()), n_(std::sqrt(inpt.size()))
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

    template <std::forward_iterator RandIt>
    matrix_t(RandIt start, RandIt fin)
        : internal::matrix_buff<T>(std::distance(start, fin)),
          n_(std::sqrt(std::distance(start, fin)))
    {
        for (std::size_t i = 0; start != fin; ++i, ++start)
        {
            this->construct(data_ + i, *start);
            ++used_;
        }
        assert(used_ == sz_);
    }

    matrix_t &operator=(const std::initializer_list<T> &inpt)
    {
        if (inpt.size() != sz_)
            throw MatrExcepts::wrong_init_list(
                "Wrong size of initializer list.");
        std::move(inpt.cbegin(), inpt.cend(), data_);
        assert(used_ == sz_);
        return *this;
    }

    iterator begin() noexcept { return iterator{data_}; }
    const_iterator cbegin() const noexcept { return const_iterator{data_}; }
    iterator end() noexcept { return iterator{data_ + sz_}; }
    const_iterator cend() const noexcept { return const_iterator{data_ + sz_}; }

    virtual ~matrix_t() = default;

public:
    matrix_t &operator*=(const T &val) noexcept
    {
        std::transform(begin(), end(), begin(),
                       [&val](auto &el) { return el * val; });
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
        auto start = iterator{access(row_id, 0)},
             fin = iterator{access(row_id + 1, 0)};
        std::transform(start, fin, start,
                       [&val](auto &el) { return el * val; });
    }

    void div_row(std::size_t row_id, const T &val) noexcept
    {
        auto start = iterator{access(row_id, 0)},
             fin = iterator{access(row_id + 1, 0)};
        std::transform(start, fin, start,
                       [&val](auto &el) { return el / val; });
    }

    void add_row(std::size_t lhs, std::size_t rhs, const T &val = T(1)) noexcept
    {
        auto start1 = iterator{access(lhs, 0)},
             fin1 = iterator{access(lhs + 1, 0)},
             start2 = iterator{access(rhs, 0)},
             fin2 = iterator{access(rhs + 1, 0)};
        std::ranges::transform(
            start1, fin1, start2, fin2, start1,
            [&val](auto &el1, auto &el2) { return el1 + el2 * val; });
    }

    void sub_row(std::size_t lhs, std::size_t rhs, const T &val = T(1)) noexcept
    {
        add_row(lhs, rhs, -val);
    }

protected:
    int make_non_zero_row(std::size_t row_id) noexcept
    {
        if (!DblCmp::are_eq(*access(row_id, row_id), T(0)))
            return 1;
        for (std::size_t i = row_id + 1; i < n_; ++i)
        {
            if (!DblCmp::are_eq(*access(i, row_id), T(0)))
            {
                swap_rows(row_id, i);
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

private:
    T calculate_det_echelon() const
    {
        matrix_t<T> tmp{*this};
        int det_coeff = tmp.make_echelon_form();
        return det_coeff * tmp.diag_multiplication();
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