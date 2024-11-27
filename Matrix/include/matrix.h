#pragma once

#include "matrix_exceptions.h"

#include <iostream>
#include <cassert>
#include <utility>
#include <algorithm>
#include <optional>
#include <cmath>
#include <concepts>
#include <type_traits>
#include <string>

namespace Matrices
{

template <typename T>
concept matrix_elem = requires(T t) {
    requires std::is_floating_point_v<T>;
    requires std::is_nothrow_constructible_v<T>;
    requires std::is_nothrow_move_constructible_v<T>;
};

#if 0
template <matrix_elem T>
class diagonal_input : public std::vector<T> {
public:
    diagonal_input(std::initializer_list<T> &inpt) : std::vector<T>(inpt) {}
};
#endif

template <matrix_elem T>
class matrix_buff {
protected:
    std::size_t sz_;
    std::size_t used_;
    T *data_;

protected:
    static void construct(T *p, const T &val) {new (p) T(val);}
    static void construct(T *p, T &&val) {new (p) T(std::move(val));}
    static void destroy(T *p) {p->~T();}

    explicit matrix_buff(std::size_t sz)
    : sz_(sz), used_(0),
    data_((sz_ == 0) ? nullptr : static_cast<T*> (::operator new (sizeof(T) * sz_)))
    {
    }

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
    {
    }

    matrix_buff& operator=(const matrix_buff &m)
    {
        if (&m == this) return *this;
        matrix_buff<T> tmp(m); 
        std::swap(data_, tmp.data_);
        std::swap(sz_, tmp.sz_);
        std::swap(used_, tmp.used_);
        return *this;
    }

    matrix_buff& operator=(matrix_buff &&m) noexcept
    {
        if (&m == this) return *this;
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
    using matrix_buff<T>::data_;
    using matrix_buff<T>::sz_;
    using matrix_buff<T>::used_;
    size_t n_;
    std::optional<T> det_;

public:
    typedef T* iterator;
    typedef const T* const_iterator;

private:
    T* access(std::size_t row, std::size_t col) const
    {
        assert(row < n_ && col < n_);
        return data_ + n_ * row + col;
    }

#ifndef NDEBUG
    void dump_internal() const
    {
        for (std::size_t i = 0; i < sz_; ++i)
        {
            std::cout << data_[i] << "\t";
            if (i % n_ == n_ - 1) std::cout << std::endl;
        }
    }

    void dump(const char* descr) const
    {
        std::cout << descr << std::endl;
        dump_internal();
    }

    void dump(const std::string &descr) const
    {
        std::cout << descr << std::endl;
        dump_internal();
    }
#endif // NDEBUG

public:
    std::size_t rank() const {return n_;}
    const T& det() const
    {
        if (!det_.has_value()) throw MatrExcepts::no_det("The determinant isn't calculated.");
        return det_.value();
    }

    explicit matrix_t(std::size_t n, const T &val = T(0))
        : matrix_buff<T>(n * n), n_(n)
    {
        for (std::size_t i = 0; i < sz_; ++i)
        {
            construct(data_ + i, val);
            ++used_;
        }
        assert(used_ == sz_);
    }

    matrix_t(const std::initializer_list<T> &inpt) : matrix_buff<T>(inpt.size()), n_(std::sqrt(inpt.size()))
    {
        std::size_t i = 0;
        for (auto start = inpt.begin(), fin = inpt.end(); start != fin; ++start, ++i)
        {
            this->construct(data_ + i, *start);
            ++used_;
        } 
        assert(used_ == sz_);
    }

    template <typename RandIt>
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

#if 0
    explicit matrix_t(const diagonal_input<T> &inpt) : matrix_buff<T>(inpt.size()), n_(inpt.size())
    {
        if (n_ == 0) return;

        auto inpt_back = std::prev(inpt.cend());
        std::size_t diag_n = 0;
        for (auto start = inpt.cbegin(); start != inpt_back; ++start, ++diag_n)
        {
            construct(access(diag_n, diag_n), *start);
            for (std::size_t i = 0; i < n_; ++i)
                construct(access(diag_n, i + 1), T(0));
        }
        construct(access(n_ - 1, n_ - 1), *inpt_back);

        assert(used_ == sz_);
    }
#endif

    matrix_t& operator=(const std::initializer_list<T> &inpt)
    {
        if (inpt.size() != sz_) throw MatrExcepts::wrong_init_list("Wrong size of initializer list.");
        std::move(inpt.cbegin(), inpt.cend(), data_)
        assert(used_ == sz_);
    }

    matrix_t& operator*(const T &val) noexcept 
    {
        for (std::size_t i = 0; i < sz_; ++i)
            data_[i] *= val;
    }

    iterator       begin() noexcept {return data_;}
    const_iterator begin() const noexcept {return data_;}
    iterator       end()   noexcept {return data_ + sz_;}
    const_iterator end()   const noexcept {return data_ + sz_;}

private:
    class row_t final{
        T* base_;
    
    public:
        row_t(T* base) noexcept : base_(base) {}

        const T* base() const {return base_;}

        T operator[](std::size_t off) const
        {
            return *(base_ + off);
        }

        T& operator[](std::size_t off)
        {
            return *(base_ + off);
        }
    };

    void swap_rows(std::size_t lhs, std::size_t rhs) noexcept
    {
        if (lhs == rhs) return;
        for (std::size_t i = 0; i < n_; ++i)
            std::swap(*access(lhs, i), *access(rhs, i));
    }

    void mul_row(std::size_t row_id, const T &val) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(row_id, i) *= val;
    }

    void add_row(std::size_t lhs, std::size_t rhs) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(lhs, i) += *access(rhs, i);
    }

    void sub_row(std::size_t lhs, std::size_t rhs) noexcept
    {
        for (std::size_t i = 0; i < n_; ++i)
            *access(lhs, i) -= *access(rhs, i);
    }

public:
    row_t operator[](std::size_t off) const
    {
        return row_t{access(off), n_};
    }

private:
    T make_non_zero_row(std::size_t row_id) noexcept
    {
        for (std::size_t i = row_id; i < n_; ++i)
        {
            if (*access(i, row_id) != T(0))
            {
                swap_rows(row_id, i);
                return T(-1);
            }
        }
        return T(0);
    }
    
    T normalize_rows(std::size_t row_id) noexcept
    {
        assert(*access(row_id, row_id) != T(0));
        T transform_coeff{1};
        for (std::size_t i = row_id + 1; i < n_; ++i)
        {
            T norm_coeff = *access(i, row_id) / *access(row_id, row_id);
            if (norm_coeff == T(0)) continue;
            mul_row(row_id, norm_coeff);
            transform_coeff /= norm_coeff;
            sub_row(i, row_id); 
        }
        return transform_coeff;
    }

    T make_echelon_form()
    {
        T transform_coeff{1};
        for (std::size_t i = 0; i < n_; ++i)
        {
            transform_coeff *= make_non_zero_row(i); 
            if (transform_coeff == T(0)) return transform_coeff;
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

public:
    const T& calculate_det()
    {
        if (det_.has_value()) return det_.value();
        matrix_t<T> tmp{*this};
        T det_coeff = tmp.make_echelon_form();
        det_ = det_coeff * tmp.diag_multiplication();
        return det_.value();
    }
};

template <matrix_elem T>
struct matrix_dumper {
    void operator()(const matrix_t<T> &m) const
    {
        auto rank = m.rank(), i = 0;
        std::cout << std::endl;
        for (auto start = m.begin(), fin = m.end(); start != fin; ++start, ++i)
        {
            std::cout << *start << "\t";
            if (i % rank == rank - 1) std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};

}