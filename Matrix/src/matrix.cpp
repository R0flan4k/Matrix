#include "matrix.h"

#include <iostream>
#include <vector>

int main()
{
    std::vector<float> v{0, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrices::matrix_t<double> matr{v.cbegin(), v.cend()};
    std::cout << matr.calculate_det() << std::endl;
}