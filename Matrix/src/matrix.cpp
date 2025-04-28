#include "matrix_wrap.h"

#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

using matrix_t = typename Matrices::matrix_wrap_t<double>;

int main()
{
    int rank;
    std::cin >> rank;
    if (!std::cin.good() || rank < 0)
    {
        std::cerr << "Matrix rank input error.\n";
        return 1;
    }
    if (rank == 0)
        return 0;

    std::vector<double> input(std::istream_iterator<double>(std::cin),
                              std::istream_iterator<double>());
    if (input.size() != rank * rank)
    {
        std::cerr << "Matrix data input size error.\n";
        return 1;
    }

    try
    {
        const matrix_t matr(input.cbegin(), input.cend());
        std::cout << matr.calculate_det() << std::endl;
    } catch (const MatrExcepts::matrix_error &me)
    {
        std::cerr << me.what() << '\n';
        return 1;
    }

    return 0;
}