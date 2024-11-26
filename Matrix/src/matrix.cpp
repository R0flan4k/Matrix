#include "matrix.h"

#include <iostream>

int main()
{
    Matrices::matrix_t<double> matr{0, 2, 3, 4, 5, 6, 7, 8, 9};
    std::cout << matr.calculate_det() << std::endl;
}