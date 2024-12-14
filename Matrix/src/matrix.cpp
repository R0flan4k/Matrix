#include "matrix.h"

#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

using Matrices::const_matrix_t;

int main()
{
    size_t rank;
    std::cin >> rank;
    if (!std::cin.good())
        throw std::runtime_error("Matrix rank input error.");
    if (rank == 0)
        return 0;

    std::vector<double> input(std::istream_iterator<double>(std::cin),
                              std::istream_iterator<double>());
    if (input.size() != rank * rank)
        throw std::runtime_error("Matrix data input size error.");

    try
    {
        const_matrix_t<double> matr(input.cbegin(), input.cend());
        std::cout << matr.calculate_det() << std::endl;
    } catch (const MatrExcepts::matrix_error &me)
    {
        std::cerr << me.what() << '\n';
        return 1;
    }

    return 0;
}