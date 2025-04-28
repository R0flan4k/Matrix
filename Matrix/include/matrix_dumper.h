#pragma once

#include "matrix_data.h"

#include <iostream>

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

} // namespace Matrices