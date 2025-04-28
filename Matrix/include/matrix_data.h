#pragma once

#include "matrix_exceptions.h"
#include "matrix_iterator.h"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace Matrices {

namespace internal {

template <typename T> concept matrix_elem = requires(T t)
{
    requires std::is_convertible_v<T, double>;
    requires std::is_nothrow_move_constructible_v<T>;
};

template <typename T> concept integral_elem = matrix_elem<T> and requires(T t)
{
    requires std::is_integral_v<T>;
};

template <class U, class T> concept t_like = std::is_convertible<U, T>::value;

template <typename T> class matrix_buff {
protected:
    std::size_t sz_;
    std::size_t used_;
    T *data_;

protected:
    template <t_like<T> U, typename K> static void construct(K *p, U &&val)
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
        try
        {
            for (std::size_t i = 0; i < sz_; ++i)
            {
                construct(data_ + i, m.data_[i]);
                ++used_;
            }
        } catch (...)
        {
            this->~matrix_buff();
            throw;
        }
    }

    matrix_buff(matrix_buff &&m) noexcept
        : sz_(m.sz_), used_(m.used_), data_(m.data_)
    {
        m.data_ = nullptr; // To avoid double global operator delete.
    }

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

template <matrix_elem T> class matrix_data : private matrix_buff<T> {
protected:
    using matrix_buff<T>::data_;
    using matrix_buff<T>::sz_;
    using matrix_buff<T>::used_;
    size_t n_;

public:
    typedef random_access_iterator<T> iterator;
    typedef random_access_iterator<const T> const_iterator;

protected:
    T *access(std::size_t row, std::size_t col) const
    {
        assert(row < n_ + 1 && col < n_ + 1);
        return data_ + n_ * row + col;
    }

public:
    std::size_t rank() const { return n_; }

    explicit matrix_data(std::size_t n, const T &val = T(0))
        : matrix_buff<T>(n * n), n_(n)
    {
        for (std::size_t i = 0; i < sz_; ++i)
        {
            this->construct(data_ + i, val);
            ++used_;
        }
        assert(used_ == sz_);
    }

    matrix_data(const std::initializer_list<T> &inpt)
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

    template <std::forward_iterator FwdIt>
    matrix_data(FwdIt start, FwdIt fin)
        : matrix_buff<T>(std::distance(start, fin)),
          n_(std::sqrt(std::distance(start, fin)))
    {
        for (std::size_t i = 0; start != fin; ++i, ++start)
        {
            this->construct(data_ + i, *start);
            ++used_;
        }
        assert(used_ == sz_);
    }

    template <typename U>
    requires std::is_convertible_v<U, T> matrix_data &
    operator=(const std::initializer_list<U> &inpt)
    {
        if (inpt.size() != sz_)
            throw MatrExcepts::wrong_init_list(
                "Wrong size of initializer list.");
        std::move(inpt.begin(), inpt.end(), data_);
        assert(used_ == sz_);
        return *this;
    }

    iterator begin() noexcept { return iterator{data_}; }
    const_iterator cbegin() const noexcept { return const_iterator{data_}; }
    iterator end() noexcept { return iterator{data_ + sz_}; }
    const_iterator cend() const noexcept { return const_iterator{data_ + sz_}; }

public:
    matrix_data &operator*=(const T &val) noexcept
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

    T trace() const noexcept
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

public:
    matrix_data(const matrix_data &) = default;

    matrix_data(matrix_data &&) noexcept = default;

    matrix_data &operator=(const matrix_data &oth)
    {
        if (&oth == this)
            return *this;
        n_ = oth.rank();
        matrix_buff<T>::operator=(oth);
        return *this;
    }

    matrix_data &operator=(matrix_data &&oth) noexcept
    {
        if (&oth == this)
            return *this;
        n_ = oth.rank();
        matrix_buff<T>::operator=(std::move(oth));
        return *this;
    }

    virtual ~matrix_data() = default;
};

} // namespace internal

} // namespace Matrices