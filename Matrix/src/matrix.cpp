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

    std::vector<float> input(std::istream_iterator<float>(std::cin),
                             std::istream_iterator<float>());
    if (input.size() != rank * rank)
        throw std::runtime_error("Matrix data input size error.");

    const_matrix_t<float> matr(input.cbegin(), input.cend());
    std::cout << matr.calculate_det() << std::endl;
}