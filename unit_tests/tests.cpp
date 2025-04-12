#define private public
#define protected public
#include "matrix.h"
#undef private
#undef protected

#include "double_comparing.h"
#include "gtest/gtest.h"

#include <initializer_list>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

using Matrices::matrix_t;

TEST(Matrices, Determinant)
{
    const matrix_t<float> m1{2}, m2{0}, m3{2, 13, -6, 4}, m4{0, 13, -6, 4},
        m5{2, 0, -6, 4}, m6{2, 13, 0, 4}, m7{2, 13, -6, 0}, m8{0, 0, -6, 4},
        m9{0, 13, 0, 4}, m10{0, 13, -6, 0}, m11{2, 0, 0, 4}, m12{2, 0, -6, 0},
        m13{2, 13, 0, 0}, m14{0, 0, 0, 4}, m15{0, 0, -6, 0}, m16{0, 13, 0, 0},
        m17{2, 0, 0, 0}, m18(2), m19{25, 0, 1, 4, -8, 1, 4, 3, 1};

    EXPECT_EQ(m1.calculate_det(), 2);
    EXPECT_EQ(m2.calculate_det(), 0);
    EXPECT_EQ(m3.calculate_det(), 86);
    EXPECT_EQ(m4.calculate_det(), 78);
    EXPECT_EQ(m5.calculate_det(), 8);
    EXPECT_EQ(m6.calculate_det(), 8);
    EXPECT_EQ(m7.calculate_det(), 78);
    EXPECT_EQ(m8.calculate_det(), 0);
    EXPECT_EQ(m9.calculate_det(), 0);
    EXPECT_EQ(m10.calculate_det(), 78);
    EXPECT_EQ(m11.calculate_det(), 8);
    EXPECT_EQ(m12.calculate_det(), 0);
    EXPECT_EQ(m13.calculate_det(), 0);
    EXPECT_EQ(m14.calculate_det(), 0);
    EXPECT_EQ(m15.calculate_det(), 0);
    EXPECT_EQ(m16.calculate_det(), 0);
    EXPECT_EQ(m17.calculate_det(), 0);
    EXPECT_EQ(m18.calculate_det(), 0);
    EXPECT_EQ(m19.calculate_det(), -231);
}

TEST(Matrices, Determinant1)
{
    size_t rank = 100;
    std::string input_str{
        "42 0 1 1 1 0 0 1 0 1 1 1 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0 1 1 1 1 0 1 0 "
        "0 0 0 0 1 0 1 0 0 1 0 0 1 0 0 1 1 0 1 0 1 1 0 1 0 1 0 1 0 0 1 0 0 0 0 "
        "0 0 1 0 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 1 0 1 0 0 1 0 1 1 0 1 1 0 1 0 0 "
        "1 1 0 1 1 1 1 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 1 0 1 1 0 0 0 1 0 1 0 1 1 "
        "1 0 0 0 1 1 0 0 0 1 0 1 1 0 0 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 1 1 1 1 "
        "0 0 0 0 0 1 1 1 0 1 0 0 0 0 1 0 1 1 1 0 1 1 0 1 1 1 42 0 2 1 2 0 0 1 "
        "1 1 1 2 2 2 0 0 0 0 1 0 0 0 2 0 1 0 0 1 2 2 1 1 2 0 1 0 0 1 1 1 1 1 1 "
        "2 0 1 2 0 1 1 1 1 1 1 1 1 1 1 0 2 1 1 0 0 1 1 0 1 1 1 0 2 0 1 0 0 1 1 "
        "2 0 2 1 1 0 1 0 0 0 1 0 2 1 0 2 0 2 2 0 1 1 0 1 0 1 2 1 1 2 1 2 1 1 1 "
        "1 1 0 1 1 0 2 1 0 1 1 0 2 2 0 2 2 1 0 1 2 0 2 0 2 1 1 0 0 1 1 2 0 0 0 "
        "1 0 2 2 0 1 2 0 1 1 1 1 1 1 1 0 0 0 1 0 2 1 1 2 1 2 1 0 1 0 0 2 1 2 1 "
        "2 1 0 1 0 1 0 1 2 2 0 2 2 0 2 2 2 42 1 2 2 5 2 2 3 3 3 3 4 4 4 2 1 1 "
        "2 1 3 1 0 3 1 1 3 2 2 4 4 2 2 3 3 1 3 0 4 3 2 2 1 3 3 2 2 3 0 2 2 4 3 "
        "2 3 3 2 2 2 1 3 2 3 2 0 1 2 2 2 3 2 2 5 1 3 1 0 3 2 2 3 4 4 3 2 3 1 1 "
        "1 3 0 3 4 3 3 3 4 3 2 4 4 0 1 0 1 3 3 2 2 3 3 2 3 3 2 2 1 2 3 1 4 1 0 "
        "1 1 0 3 2 1 3 2 2 2 1 3 1 4 0 3 3 2 2 1 2 1 3 1 1 1 1 1 3 2 1 2 2 1 1 "
        "1 1 1 1 2 3 0 0 1 3 1 2 1 3 4 1 3 1 1 2 1 0 3 2 3 2 3 3 1 2 2 2 0 2 4 "
        "3 1 4 2 2 3 4 4 42 0 1 2 3 2 3 3 2 4 3 4 4 4 2 2 2 3 2 3 1 1 2 1 0 2 "
        "1 2 3 3 3 2 3 3 2 3 0 3 4 1 3 2 2 2 2 1 2 2 0 2 4 1 3 2 2 2 2 2 1 2 1 "
        "3 3 1 1 2 4 1 1 1 2 4 0 3 2 1 3 2 2 2 2 3 3 3 4 2 3 2 2 0 3 3 3 3 3 3 "
        "3 3 5 4 0 1 1 1 4 3 3 4 4 5 3 4 5 5 2 2 3 3 1 4 2 1 3 1 2 4 3 2 5 4 2 "
        "3 4 4 4 5 1 5 5 3 2 4 3 3 3 3 2 3 3 1 5 3 3 4 2 1 4 3 1 3 2 3 5 2 0 3 "
        "5 3 4 3 4 6 1 5 3 1 5 3 3 3 3 6 2 5 5 2 3 3 2 0 4 5 4 4 5 4 4 4 6 4 "
        "42 0 2 2 3 1 2 4 3 6 3 5 5 6 2 2 2 3 3 2 3 1 3 2 2 2 3 2 5 5 3 2 6 3 "
        "4 4 1 3 4 3 3 5 2 3 2 3 2 3 2 2 4 2 4 3 3 2 5 3 2 4 3 3 3 2 1 2 5 3 3 "
        "3 2 5 1 4 3 2 4 2 4 1 3 4 2 4 5 2 3 2 2 1 5 3 3 5 3 5 5 4 5 3 0 1 0 0 "
        "1 1 0 0 0 -1 0 -2 -1 -4 -3 -4 -1 -3 -4 -3 -6 -6 -5 -8 -9 -6 -8 -8 -9 "
        "-5 -6 -7 -10 -9 -7 -12 -11 -5 -8 -12 -13 -9 -4 -12 -7 -11 -8 -8 -11 "
        "-11 -9 -10 -12 -8 -9 -9 -14 -8 -12 -11 -11 -9 -12 -8 -9 -13 -9 -10 -5 "
        "-9 -7 -11 -4 -7 -9 -10 -13 -9 -7 -13 -9 -10 -9 -10 -15 -9 -11 -11 -5 "
        "-13 -9 -8 -7 -10 -9 -9 -14 -11 -8 -12 42 1 1 2 3 2 1 4 3 6 4 5 5 5 2 "
        "2 2 3 3 3 3 1 2 2 2 4 5 3 5 4 3 2 5 4 2 5 1 3 5 3 3 3 2 3 4 3 1 3 2 2 "
        "4 2 4 2 4 2 4 4 2 4 2 3 4 2 1 1 4 2 3 3 3 5 3 4 3 3 2 1 2 2 4 3 3 4 3 "
        "2 3 3 4 2 5 4 4 3 4 5 5 4 6 4 0 0 1 1 3 2 2 1 4 4 3 6 4 5 4 3 2 5 3 5 "
        "4 1 2 3 3 3 4 2 4 4 3 4 5 4 2 5 1 4 3 4 5 3 3 3 3 4 2 1 1 4 4 3 3 4 4 "
        "3 4 1 4 4 4 3 2 1 2 2 6 3 3 3 2 6 2 3 1 3 4 3 2 4 4 4 3 3 7 2 4 3 3 3 "
        "5 6 4 5 3 4 6 4 5 4 42 0 2 1 3 2 2 2 4 4 5 7 8 8 3 6 3 4 4 4 2 3 4 3 "
        "3 3 3 3 3 4 3 5 6 5 4 4 1 4 5 3 6 5 3 5 4 5 4 3 2 5 5 2 6 2 3 4 6 3 2 "
        "7 3 5 5 2 3 4 5 3 1 4 3 6 2 5 3 5 3 6 5 3 5 4 4 3 7 4 5 3 4 3 8 5 4 6 "
        "3 6 7 5 6 4 0 1 1 0 3 2 2 3 3 4 4 2 5 6 3 4 5 3 2 4 3 2 4 1 3 6 4 4 5 "
        "4 1 2 4 5 4 4 3 6 6 4 4 5 3 5 3 6 4 2 5 4 7 4 5 4 2 2 4 4 1 6 1 4 6 4 "
        "3 3 4 3 5 4 5 6 3 4 4 2 6 4 5 4 4 7 1 5 6 3 4 4 4 2 5 5 5 6 4 3 5 4 5 "
        "2 0 0 1 1 2 1 1 1 2 2 1 5 3 4 2 3 1 2 1 3 2 3 2 2 4 2 3 2 4 3 3 4 4 3 "
        "3 3 1 3 3 3 3 3 3 3 4 3 2 2 2 2 2 3 2 2 2 1 5 1 3 3 4 1 2 1 2 3 3 2 3 "
        "3 1 4 1 4 2 3 3 4 2 3 3 3 3 2 5 1 4 2 2 3 4 5 3 2 2 5 5 3 4 4 42 0 2 "
        "1 2 0 1 2 1 2 3 4 4 6 0 5 3 0 2 1 2 4 5 1 4 3 4 3 3 4 3 4 5 2 2 3 0 4 "
        "5 3 4 3 3 6 3 3 4 2 2 4 4 2 4 1 3 3 5 3 2 7 3 5 4 4 3 4 2 1 3 3 2 3 2 "
        "3 4 3 3 4 5 2 4 2 3 1 5 4 5 2 4 3 5 4 4 6 0 6 5 3 5 3 0 0 0 1 2 2 3 3 "
        "2 5 4 4 5 7 4 6 7 5 5 5 4 4 3 3 5 9 7 5 6 4 3 5 6 7 4 7 3 6 9 5 6 5 5 "
        "5 5 7 5 4 3 5 7 4 8 5 4 2 6 5 4 7 3 5 9 6 4 5 7 3 5 4 6 6 4 6 8 4 7 4 "
        "4 6 4 6 4 7 9 6 8 7 5 4 6 7 7 6 6 5 7 5 9 5 0 0 0 0 1 2 2 2 3 5 5 4 5 "
        "5 3 4 3 5 4 4 4 3 2 2 4 4 5 5 4 3 2 3 6 5 4 5 3 5 6 3 4 5 2 4 2 5 1 5 "
        "2 4 5 2 6 3 2 3 4 5 2 5 2 4 6 5 2 4 6 4 2 5 4 5 3 3 3 3 4 3 2 3 4 5 2 "
        "4 6 4 4 6 4 3 5 5 5 6 3 4 6 4 6 3 0 1 1 0 2 1 0 1 2 1 2 3 4 4 2 5 3 2 "
        "4 4 2 4 4 3 6 4 5 4 4 5 3 5 4 5 4 4 4 5 4 4 5 4 2 6 4 5 5 3 5 5 5 6 4 "
        "2 3 3 6 3 3 7 6 5 5 4 4 5 1 3 4 5 4 6 2 4 3 4 3 5 3 7 5 5 3 3 7 3 4 3 "
        "5 3 5 6 5 4 2 8 4 4 5 5 0 1 0 0 1 1 1 3 1 4 4 1 4 4 2 3 4 2 3 4 3 3 2 "
        "2 4 5 4 5 3 5 1 3 5 6 3 2 5 5 5 3 3 3 1 6 3 5 2 4 4 2 7 3 6 4 3 2 5 5 "
        "1 5 2 3 6 5 2 3 3 1 4 4 4 6 3 5 4 3 5 3 4 4 3 5 2 6 5 3 5 4 4 3 5 5 6 "
        "4 4 4 3 4 6 2 0 0 1 0 2 1 2 2 3 4 4 3 5 7 4 6 5 4 4 4 6 3 4 3 5 6 5 6 "
        "5 7 1 4 7 5 4 6 5 5 6 7 6 7 3 6 3 8 5 3 5 6 8 5 8 5 4 4 7 3 4 8 5 6 8 "
        "6 4 5 5 4 6 5 6 7 3 4 5 4 8 5 6 5 6 7 2 5 10 6 7 5 6 5 6 6 7 9 4 4 6 "
        "5 5 4 0 1 0 1 2 1 2 4 1 5 3 1 2 4 3 4 6 2 4 3 5 5 6 2 5 7 7 6 9 8 4 3 "
        "7 7 4 7 6 7 9 8 5 6 3 8 4 7 5 4 4 4 8 7 7 5 6 3 5 5 6 8 6 8 8 8 5 5 5 "
        "4 7 7 4 6 4 4 6 3 10 3 3 9 5 8 4 8 10 6 7 7 7 6 5 7 9 5 5 8 6 5 8 7 "
        "42 0 2 1 2 0 0 2 2 4 4 4 4 5 3 3 2 2 5 2 4 5 5 4 8 5 6 8 8 7 3 6 9 6 "
        "4 4 6 5 7 7 6 6 3 9 1 8 5 3 4 5 6 5 10 6 6 5 7 5 4 8 6 5 5 5 4 6 4 7 "
        "5 7 3 6 4 4 4 3 8 4 4 5 6 6 4 5 9 4 4 6 6 7 7 6 7 6 5 9 9 4 6 6 42 1 "
        "1 1 2 2 1 3 2 4 3 3 5 4 1 4 6 3 7 6 4 6 9 7 7 6 7 7 9 7 7 7 9 6 5 9 6 "
        "5 7 8 8 7 1 12 5 5 5 6 5 7 10 7 8 3 8 7 7 5 5 12 7 10 11 7 6 7 4 4 3 "
        "6 7 10 3 6 5 6 7 5 5 10 5 7 6 9 12 5 9 6 5 5 10 8 6 6 8 10 6 10 9 12 "
        "0 1 1 0 2 1 0 1 3 3 3 3 4 5 4 6 5 4 6 6 6 5 6 6 9 7 9 7 7 9 5 6 6 7 4 "
        "8 7 7 5 8 10 7 3 10 4 8 6 2 7 9 10 8 7 4 7 7 8 3 5 10 9 8 6 5 6 6 4 6 "
        "7 6 7 9 4 5 5 6 7 6 6 9 6 8 3 5 13 5 7 5 8 7 7 9 8 8 5 10 7 8 5 8 42 "
        "0 2 1 2 0 0 1 2 3 4 4 5 6 3 4 4 3 6 2 3 3 4 3 4 6 5 5 6 6 3 4 6 7 3 4 "
        "3 6 5 6 7 4 4 8 4 8 6 3 4 6 5 3 6 5 4 4 5 5 4 8 4 6 4 4 5 4 4 3 4 5 3 "
        "5 5 4 5 6 5 3 5 5 6 4 5 5 7 6 5 6 6 5 8 6 6 6 2 7 8 3 6 2 42 1 1 2 4 "
        "2 3 5 3 7 7 6 7 9 6 7 5 6 7 7 7 7 9 8 10 11 13 11 13 10 8 7 13 14 7 "
        "12 9 10 11 9 12 9 5 14 7 13 7 7 7 11 15 9 12 7 11 9 10 9 8 15 9 12 11 "
        "9 9 10 9 9 6 11 7 13 7 8 9 8 13 9 6 14 10 13 10 12 16 8 11 11 11 10 "
        "13 12 11 9 11 15 13 11 13 13 42 1 1 2 4 3 3 4 4 6 5 5 7 6 5 6 6 6 8 6 "
        "4 4 7 7 6 8 8 9 10 9 7 7 8 10 5 12 5 6 9 9 10 9 4 10 6 8 8 7 5 9 10 9 "
        "10 4 7 6 9 6 8 11 8 10 11 6 6 8 7 6 4 7 8 10 5 6 6 7 9 7 3 12 9 9 7 "
        "10 13 8 8 6 8 4 10 8 8 7 8 11 8 10 10 13 0 0 1 0 2 2 1 1 3 3 2 3 3 3 "
        "3 2 3 4 4 4 3 4 5 4 8 7 6 8 12 5 6 7 8 6 8 9 6 6 8 9 8 8 4 9 4 8 6 7 "
        "6 8 6 7 9 6 5 5 8 6 7 7 4 5 8 6 4 10 8 9 5 7 8 8 5 3 5 4 10 5 4 10 8 "
        "9 4 7 10 5 7 8 4 4 8 7 6 8 7 8 11 7 8 9 0 1 0 0 2 3 1 1 3 2 3 4 5 5 3 "
        "6 5 3 3 6 4 6 6 5 7 9 8 7 9 7 4 8 5 9 5 11 5 6 8 9 11 7 5 11 9 10 6 5 "
        "7 8 10 6 8 6 6 5 7 5 6 10 7 8 9 4 7 8 6 5 6 6 9 8 5 9 9 8 8 8 7 9 7 9 "
        "8 9 12 8 11 10 8 10 10 12 9 6 9 10 11 8 9 11 42 1 1 1 3 2 2 4 2 4 5 3 "
        "6 6 2 5 4 2 3 3 1 3 4 2 5 8 6 5 6 6 6 5 7 8 6 7 6 7 9 5 7 4 4 7 5 7 5 "
        "6 4 7 9 5 9 2 7 6 7 6 4 8 2 8 10 7 2 8 6 6 4 6 6 8 5 6 7 4 7 5 7 8 7 "
        "7 6 5 7 7 8 5 7 3 8 5 8 6 8 9 8 6 10 7 0 0 0 1 2 2 3 2 2 3 2 3 3 3 2 "
        "2 3 3 2 4 3 4 5 5 4 5 5 5 8 4 5 8 8 8 5 8 5 5 8 9 6 5 5 8 6 5 4 6 4 6 "
        "8 6 7 6 7 4 6 4 5 7 5 7 10 4 3 8 7 6 4 6 5 7 4 8 5 4 10 5 4 6 4 7 7 8 "
        "10 6 8 7 4 5 8 7 8 5 11 7 7 8 9 9 0 0 1 1 3 1 3 2 3 4 4 5 5 7 6 6 4 5 "
        "4 6 6 6 6 6 7 6 7 5 8 10 5 7 11 11 7 10 8 9 9 10 10 8 7 11 7 11 7 5 5 "
        "8 12 10 11 8 10 6 9 5 10 12 11 10 9 6 8 10 10 9 7 10 5 10 6 10 8 5 12 "
        "8 9 10 7 11 8 9 15 9 12 8 7 11 11 12 12 8 9 12 13 9 9 10 0 1 0 0 2 3 "
        "1 2 2 3 3 4 4 5 3 5 5 2 1 7 5 4 5 4 9 9 9 6 10 7 1 8 6 8 4 9 4 1 4 5 "
        "8 3 0 5 10 10 1 -1 3 6 8 1 4 0 0 -1 0 -1 -1 0 -3 -3 1 -3 4 -3 0 -4 -2 "
        "-6 -2 -1 -4 -1 -2 -3 -7 -1 -5 1 -5 -3 -4 -1 -2 -9 0 -8 -7 -2 -1 3 -14 "
        "-5 -6 -8 -4 -8 -2 -5 42 1 1 1 2 1 0 2 1 3 3 1 3 3 2 3 5 3 5 4 3 3 6 3 "
        "4 7 7 5 9 5 4 4 6 7 6 9 7 9 10 8 12 5 4 11 8 5 7 6 7 11 9 11 9 4 9 4 "
        "8 7 5 12 4 9 9 6 8 8 5 5 6 7 9 9 5 6 8 5 6 3 7 8 7 9 5 7 9 7 8 8 8 5 "
        "8 7 10 8 6 8 9 8 7 7 42 0 2 2 3 1 2 3 2 4 3 5 6 7 2 5 5 3 5 4 4 5 6 3 "
        "4 6 4 6 8 9 7 7 8 7 8 9 8 10 9 10 9 7 6 11 8 8 6 8 6 6 9 9 8 6 6 5 10 "
        "5 6 10 6 10 10 8 6 7 5 4 7 7 6 9 3 9 7 7 8 5 9 8 6 9 8 8 11 9 12 9 7 "
        "6 9 9 11 7 6 10 10 7 10 9 0 0 -1 -1 -2 -1 0 -1 -1 0 2 -3 0 0 2 3 3 2 "
        "2 0 2 1 3 1 -2 0 1 0 1 1 1 -2 2 3 2 1 3 0 1 0 2 0 0 1 2 -2 -4 -2 -4 7 "
        "1 1 -2 -7 -4 1 -4 -3 -3 0 -4 0 0 -3 -2 -6 -3 -5 -10 -5 -8 -5 -1 -6 -8 "
        "-3 -10 -5 -5 -8 -9 -5 -5 -10 -8 -4 -11 -10 -7 -10 -5 -8 -5 -4 -9 -8 "
        "-10 -6 -12 -11 0 0 0 0 1 2 2 1 2 3 3 2 3 3 2 3 4 4 3 4 3 4 4 2 4 5 5 "
        "6 6 4 8 4 5 6 5 10 7 10 10 10 9 7 6 9 4 8 4 8 5 7 8 7 8 5 6 8 5 6 8 9 "
        "5 9 11 10 4 7 7 6 5 7 6 7 6 6 4 6 9 6 5 10 5 9 7 7 9 9 11 11 5 7 6 8 "
        "8 7 7 6 12 9 9 11 42 0 2 1 1 -1 -1 2 1 3 1 3 3 4 2 2 1 1 4 2 4 3 5 6 "
        "7 3 5 4 8 6 0 5 9 6 5 4 5 -1 3 6 7 5 -1 8 5 5 6 2 8 6 5 6 7 5 6 2 11 "
        "6 5 7 7 3 4 0 7 8 5 5 3 6 3 9 0 5 6 5 7 5 5 5 8 5 3 7 10 1 3 2 3 7 9 "
        "4 2 7 4 7 6 5 4 4 0 1 0 0 1 2 1 2 3 4 4 4 6 5 3 6 5 3 5 4 3 6 7 6 7 8 "
        "8 7 9 8 7 8 10 12 7 13 9 10 12 13 14 10 8 16 10 11 8 13 7 10 10 10 13 "
        "5 12 7 14 7 11 14 10 13 13 8 8 11 10 9 5 11 9 11 8 11 8 13 11 12 9 13 "
        "9 13 11 11 13 13 14 11 11 8 14 10 12 6 10 13 14 12 10 13 0 0 1 0 2 2 "
        "2 2 4 5 5 6 6 7 3 6 5 4 6 4 4 5 5 4 8 7 9 6 9 8 6 10 11 10 7 11 8 10 "
        "13 12 16 8 5 11 7 15 10 11 7 10 11 8 15 11 11 10 13 9 14 15 10 11 13 "
        "14 7 13 10 9 7 13 8 12 8 11 11 9 15 7 7 14 12 9 11 12 16 13 13 13 9 "
        "13 13 13 13 12 8 13 16 11 16 11 0 0 1 0 0 -2 -1 0 -1 -2 -2 -1 0 1 -1 "
        "1 1 -2 1 -1 0 1 2 3 2 1 0 -1 0 2 -4 3 1 -1 1 1 -2 -2 -1 0 1 3 0 3 3 0 "
        "7 1 4 0 3 5 4 2 0 -3 6 3 3 3 4 1 2 -1 4 5 0 2 2 1 4 1 0 3 8 -1 0 0 6 "
        "0 2 0 -1 3 3 3 0 -3 -3 2 4 1 2 3 -1 5 -1 1 0 -1 0 1 0 1 3 3 2 2 4 5 4 "
        "4 4 4 4 4 6 6 7 6 5 4 6 5 6 9 11 6 10 8 8 8 6 11 6 16 10 14 14 13 17 "
        "8 9 12 9 13 10 10 8 12 14 13 13 11 13 10 10 10 14 15 11 14 14 12 10 "
        "11 9 11 9 14 10 13 8 14 12 9 13 5 7 16 9 13 13 13 18 14 12 15 11 13 "
        "11 15 15 9 11 15 16 12 13 15 0 1 1 0 3 2 2 2 3 3 4 3 4 5 3 6 5 3 3 6 "
        "6 7 9 4 8 6 7 8 9 9 4 7 9 7 6 12 11 12 12 14 13 10 6 15 6 11 9 9 11 "
        "12 11 15 14 9 12 8 12 7 12 16 12 12 13 12 9 11 7 9 11 14 10 12 4 11 6 "
        "6 13 12 9 14 11 16 9 8 17 12 12 12 11 11 9 13 13 12 7 12 13 11 11 16 "
        "0 1 1 1 4 3 2 3 4 5 3 5 5 6 5 5 5 5 5 8 5 3 4 5 9 9 8 7 8 9 5 9 8 9 5 "
        "11 10 9 10 11 15 8 7 11 9 12 12 8 11 10 14 13 14 11 13 7 14 4 12 11 "
        "13 12 10 9 9 10 11 11 11 7 13 15 5 13 9 8 13 11 10 15 8 15 8 10 16 9 "
        "13 9 10 8 10 13 11 10 12 13 11 10 13 13 0 1 0 0 2 2 1 1 2 2 3 1 3 2 3 "
        "3 4 4 4 6 4 3 4 3 5 6 6 6 6 6 4 5 6 6 5 10 9 9 11 11 10 7 6 12 8 12 "
        "13 9 12 9 12 16 12 9 16 8 10 11 13 15 10 10 15 13 13 14 8 9 8 13 12 "
        "16 8 12 10 5 12 8 9 17 13 14 6 11 18 11 15 12 10 12 10 15 11 10 9 12 "
        "13 12 10 15 42 0 2 1 3 1 2 2 3 4 6 6 7 8 4 7 4 4 5 5 5 5 7 5 7 6 6 8 "
        "7 10 5 10 12 8 6 11 10 9 12 12 15 10 6 16 9 14 14 11 11 12 12 16 18 "
        "10 15 10 18 10 16 19 13 12 13 14 11 14 11 9 7 15 9 16 7 14 10 8 14 14 "
        "11 14 16 15 10 10 19 14 16 11 12 14 13 14 13 14 8 15 15 13 15 14 0 1 "
        "1 1 3 1 1 2 3 4 3 4 3 4 4 4 4 4 5 5 5 4 6 7 10 8 12 5 11 9 6 7 9 10 6 "
        "12 7 9 8 11 12 7 5 11 7 14 11 9 10 13 11 14 11 8 18 9 12 10 16 16 13 "
        "8 10 9 10 14 12 11 10 15 10 15 10 13 12 5 13 8 8 18 15 12 9 9 21 12 "
        "12 9 11 12 12 16 10 12 10 18 14 13 14 15 0 0 0 0 0 1 0 1 1 2 1 2 3 4 "
        "2 5 5 2 4 2 2 4 2 3 6 8 6 6 7 4 5 6 4 7 7 11 6 7 12 12 11 10 7 8 8 13 "
        "10 14 9 11 8 9 13 6 10 4 10 10 13 12 6 7 14 9 8 10 8 7 9 9 11 9 9 12 "
        "12 10 10 7 11 12 11 10 11 9 11 14 12 12 8 10 10 9 11 8 9 12 15 9 12 "
        "12 42 0 1 1 1 0 1 3 0 4 3 1 3 4 3 4 4 2 5 3 4 4 4 2 5 4 3 8 7 9 6 7 "
        "10 8 8 9 15 9 13 12 13 9 6 12 9 12 13 16 10 14 16 18 21 11 17 8 15 11 "
        "16 15 12 13 16 14 11 14 15 11 9 14 13 17 9 15 13 8 16 10 14 17 14 15 "
        "11 10 19 16 17 11 13 11 13 13 16 13 12 15 14 12 15 14 42 0 1 2 2 0 2 "
        "4 0 4 2 3 3 7 2 4 3 3 3 5 6 6 7 7 7 7 8 7 8 6 9 5 11 8 6 13 11 10 9 8 "
        "13 9 8 13 9 8 7 14 9 13 17 13 17 11 13 10 14 9 13 13 11 13 18 13 9 13 "
        "12 13 6 9 10 15 7 14 15 10 14 8 16 15 6 14 13 11 16 12 17 11 6 10 11 "
        "10 12 11 14 14 14 14 12 14 0 0 0 1 2 2 3 3 2 5 3 4 3 6 4 6 6 4 5 6 6 "
        "7 5 5 9 10 9 8 10 10 8 7 10 10 7 15 9 11 12 11 16 9 8 15 8 13 11 13 "
        "10 11 17 13 20 16 16 10 13 9 17 15 14 15 16 15 10 16 17 12 9 12 12 17 "
        "7 15 18 11 18 8 16 20 12 15 17 14 22 16 20 15 10 14 14 18 19 14 15 19 "
        "19 13 19 18 42 1 2 2 4 1 2 4 2 4 4 4 5 6 2 5 5 2 5 4 4 7 8 5 9 9 9 8 "
        "10 11 7 8 9 11 7 14 10 15 12 11 14 8 10 16 9 10 11 12 10 12 18 15 18 "
        "11 15 7 15 12 14 17 15 18 14 11 9 15 13 12 12 16 13 15 8 17 16 10 16 "
        "10 15 16 13 14 15 11 18 17 15 12 12 11 15 15 18 12 12 21 16 11 15 16 "
        "0 0 0 1 1 0 2 2 0 3 2 2 1 4 3 5 4 3 4 5 6 6 4 4 6 4 6 7 6 8 8 8 8 8 5 "
        "8 9 8 9 8 11 7 5 11 10 8 12 12 8 11 14 14 12 10 11 9 13 9 15 11 12 11 "
        "14 12 9 15 14 8 9 11 11 14 8 12 13 7 16 7 11 17 10 8 9 7 18 15 12 7 8 "
        "11 9 13 13 15 7 14 9 13 13 10 0 0 0 0 1 2 2 2 3 5 5 4 6 6 4 6 6 6 7 5 "
        "5 5 6 5 7 8 8 7 8 5 3 6 9 7 6 10 8 7 10 9 10 10 5 8 6 11 6 11 5 10 12 "
        "10 18 12 9 8 12 10 14 15 9 13 18 15 10 12 10 14 6 13 13 11 7 12 14 9 "
        "13 10 9 17 10 12 12 12 17 14 14 13 10 11 13 10 11 12 11 13 11 13 12 "
        "13 42 0 2 2 3 0 2 3 1 3 2 4 4 7 3 5 3 2 2 4 4 5 5 4 5 5 5 4 7 9 6 5 7 "
        "9 8 11 9 12 11 11 16 9 8 12 10 13 12 12 12 13 17 13 15 11 14 9 15 11 "
        "17 16 12 12 14 11 14 13 12 10 11 13 10 16 8 18 17 10 16 11 18 14 14 "
        "15 14 12 18 15 16 12 10 15 13 15 16 13 10 18 18 14 13 13 42 1 2 1 4 2 "
        "2 4 4 5 6 5 8 8 4 7 5 3 6 4 4 6 7 5 8 8 8 9 11 12 6 10 12 12 9 11 11 "
        "9 14 14 14 9 7 13 8 15 13 13 13 13 17 11 20 14 16 12 19 14 18 20 16 "
        "16 17 16 11 18 9 14 11 18 12 16 11 18 13 9 19 10 11 15 16 16 15 16 19 "
        "15 18 14 15 14 17 15 19 11 14 18 17 15 20 15 0 0 0 -1 -1 0 -2 -1 0 0 "
        "1 1 1 2 2 4 2 1 4 2 3 2 2 2 8 7 6 7 5 6 2 6 4 7 4 10 11 6 8 8 15 5 4 "
        "12 7 13 10 11 12 12 9 13 13 12 12 5 15 11 16 17 13 10 11 14 11 12 8 9 "
        "8 13 9 18 4 12 13 12 11 12 12 15 13 13 13 15 14 15 12 16 12 15 11 12 "
        "10 9 6 12 15 8 13 11 0 0 0 0 1 2 1 0 2 2 2 2 2 2 2 2 3 3 2 4 4 2 3 3 "
        "4 6 4 8 7 4 4 5 5 4 2 9 3 5 6 8 8 5 4 9 7 7 6 8 7 9 11 8 9 8 10 6 9 7 "
        "12 12 8 10 10 10 8 10 11 6 5 6 11 12 9 10 10 8 12 8 8 9 12 9 9 10 13 "
        "9 14 9 8 9 10 11 11 10 10 8 9 12 9 10 0 0 0 0 0 1 0 1 2 3 2 4 4 6 2 6 "
        "4 4 5 5 7 6 5 6 8 7 10 8 9 7 10 10 12 9 9 16 13 11 15 14 17 10 9 16 "
        "13 13 11 15 12 15 17 16 18 12 17 13 17 12 19 25 18 16 22 19 9 18 18 "
        "11 11 20 16 22 13 19 19 16 18 14 19 17 17 12 18 13 23 21 20 18 13 16 "
        "20 20 22 17 16 17 20 21 15 18 0 1 0 0 1 1 0 1 1 0 1 1 1 -1 -2 0 -1 -2 "
        "-1 -2 -3 0 2 -1 2 1 4 -1 0 0 0 3 -1 2 0 2 -2 -1 3 -3 1 1 2 4 1 -2 1 1 "
        "0 -1 1 3 0 -1 4 -2 2 3 -2 5 4 3 -1 -3 0 7 -2 -1 -5 2 -5 -1 -5 10 1 -1 "
        "-3 -2 -5 -7 -4 1 0 -1 -9 -5 -2 -2 -6 -1 -5 -9 -4 -8 -3 -2 -8 -4 -4 -6 "
        "0 1 0 0 1 2 0 1 2 3 2 2 2 1 1 1 3 3 4 5 3 3 4 2 4 4 5 4 6 6 8 5 4 7 5 "
        "10 8 11 8 10 9 7 4 10 8 7 7 13 7 8 12 10 11 9 13 10 11 15 19 16 16 17 "
        "13 13 8 14 16 11 12 17 13 16 12 16 12 11 15 11 14 17 13 14 15 13 18 "
        "16 19 17 10 13 12 13 16 10 12 15 17 13 16 15 42 0 1 1 1 1 1 3 2 5 3 3 "
        "4 4 2 4 4 3 7 3 4 6 6 5 7 5 7 8 11 9 8 10 9 8 10 11 10 7 12 13 12 9 5 "
        "11 8 9 11 11 11 14 13 11 17 12 14 10 16 12 15 18 14 15 16 14 10 17 17 "
        "13 13 16 16 17 11 14 15 9 14 11 14 15 14 15 10 15 19 14 17 14 12 12 "
        "19 14 16 15 15 23 17 17 18 17 42 0 1 1 2 2 1 2 3 4 3 4 5 6 4 6 5 4 6 "
        "4 4 5 4 4 5 7 6 6 9 7 4 6 7 7 7 10 6 6 13 12 14 10 6 9 7 10 14 10 12 "
        "13 14 9 20 13 15 10 15 14 19 21 15 15 18 16 13 19 14 13 12 17 15 18 "
        "14 16 19 11 19 11 19 13 18 18 12 15 21 16 20 20 15 16 17 16 21 13 16 "
        "20 21 16 16 16 42 0 2 2 3 1 2 4 2 5 2 4 5 7 3 5 5 4 5 5 5 5 6 6 9 9 9 "
        "9 13 11 9 9 11 9 10 14 10 10 13 15 17 13 9 15 11 14 16 16 16 16 16 16 "
        "21 16 18 10 20 15 21 22 17 15 20 19 12 22 17 17 17 19 20 22 12 19 22 "
        "12 22 15 24 21 18 19 15 18 24 20 23 19 13 13 19 17 19 20 15 22 21 21 "
        "19 21 0 1 1 1 3 1 2 3 3 4 4 5 6 8 4 7 4 3 3 5 6 7 8 7 8 9 10 9 12 13 "
        "10 11 12 15 9 15 13 14 13 17 18 12 10 19 15 15 15 16 16 15 19 17 18 "
        "16 18 14 24 14 22 22 19 17 21 19 13 23 19 16 17 23 15 21 15 21 20 17 "
        "26 19 17 23 22 21 19 18 25 21 24 18 15 19 21 19 23 17 17 25 23 23 21 "
        "20 0 1 0 0 2 2 1 1 3 3 4 3 4 5 5 7 5 5 6 7 7 5 6 7 9 9 11 11 9 9 7 9 "
        "9 11 5 14 10 10 12 14 16 10 7 16 9 15 15 12 14 16 17 15 15 10 21 11 "
        "16 13 18 23 18 17 18 14 17 21 17 15 18 20 20 25 14 18 15 13 23 18 18 "
        "19 22 19 15 17 24 21 21 19 22 20 19 19 20 18 15 19 18 23 20 22 42 1 2 "
        "2 4 1 1 3 3 4 4 5 6 6 4 5 3 3 4 3 2 2 4 5 7 8 8 5 7 9 5 6 7 10 5 10 6 "
        "8 9 9 11 7 7 11 6 11 9 8 9 11 9 15 13 8 15 6 14 10 13 19 13 11 14 11 "
        "8 14 12 13 15 17 13 20 11 16 15 9 17 14 17 14 19 16 15 13 18 17 13 14 "
        "17 15 18 15 16 11 14 20 17 14 18 21 0 0 1 0 1 1 1 2 3 4 3 5 5 7 1 6 4 "
        "3 4 3 6 5 7 5 8 6 9 8 8 6 7 9 11 6 7 16 8 9 12 12 14 13 9 15 10 11 9 "
        "17 9 11 15 14 18 9 13 9 19 10 19 21 18 17 19 16 12 17 18 18 15 18 19 "
        "18 13 18 17 14 22 15 17 14 18 16 16 15 21 21 21 18 15 14 19 17 20 16 "
        "15 19 20 18 21 21 0 1 0 1 2 1 1 2 1 3 3 3 3 5 3 4 5 4 5 6 6 5 6 5 6 8 "
        "9 7 7 7 8 7 8 11 3 12 10 12 11 11 11 7 9 15 12 12 11 12 10 11 16 16 "
        "15 14 14 12 13 14 20 20 18 17 20 16 12 11 19 11 16 22 15 23 14 19 20 "
        "13 21 18 15 22 19 15 20 16 25 21 20 16 19 17 18 23 17 14 16 20 18 17 "
        "20 18 0 0 1 0 2 1 1 1 2 2 3 3 3 5 3 4 3 2 2 3 5 5 5 3 8 7 7 9 10 6 3 "
        "10 10 8 7 9 9 9 14 12 11 8 10 13 9 14 16 13 14 13 17 13 18 14 17 8 15 "
        "15 18 21 16 11 16 15 14 21 18 18 17 25 19 23 15 19 22 10 25 17 18 13 "
        "22 19 14 16 26 20 20 22 18 17 19 22 22 17 17 19 23 15 19 15 0 1 1 1 3 "
        "2 2 4 3 5 2 4 4 5 2 4 5 3 4 6 5 6 7 6 10 7 10 5 11 9 7 11 10 9 9 11 8 "
        "8 11 11 9 7 7 10 9 8 10 9 12 11 15 10 13 11 15 8 14 13 14 15 17 12 15 "
        "11 11 18 14 17 14 16 14 20 9 18 16 11 18 13 18 15 15 16 13 14 21 13 "
        "17 14 14 16 16 16 18 13 17 23 13 16 19 18 0 0 0 0 1 1 1 1 1 1 1 2 1 3 "
        "2 3 1 1 2 2 3 3 4 4 6 6 6 7 9 6 5 9 8 10 6 12 9 7 10 10 15 7 7 11 10 "
        "11 13 17 10 14 16 11 17 15 12 10 15 13 22 16 15 13 18 15 12 21 19 17 "
        "13 21 13 21 16 16 23 14 24 15 19 19 21 20 20 17 23 19 19 18 15 16 18 "
        "18 22 12 19 23 19 17 21 16 42 1 2 2 5 2 2 4 3 4 4 6 5 8 2 4 3 3 3 6 6 "
        "7 9 7 9 9 10 10 11 6 8 9 12 8 5 14 7 11 13 10 12 9 10 17 10 12 11 13 "
        "12 12 14 14 18 13 15 13 16 15 17 21 13 11 21 13 11 18 12 13 9 14 12 "
        "24 10 18 17 12 18 11 17 17 16 17 19 15 18 19 14 17 13 17 13 20 15 15 "
        "14 19 19 18 18 18 0 1 0 1 2 2 1 2 2 4 2 2 2 2 2 2 5 3 4 4 4 2 4 1 2 5 "
        "5 4 8 8 6 5 7 6 3 10 8 9 10 10 12 6 4 10 10 9 11 11 8 6 13 14 12 10 "
        "12 8 10 7 15 13 15 12 15 15 10 14 15 11 11 16 11 19 11 15 16 9 18 6 "
        "10 17 15 17 13 14 21 15 18 18 12 16 13 18 18 10 15 16 15 14 18 14 42 "
        "1 2 1 3 1 1 3 2 4 4 2 5 4 2 3 5 3 6 4 2 3 7 4 7 8 7 7 9 10 7 6 8 9 7 "
        "11 9 10 10 10 15 10 5 17 7 11 12 13 12 11 16 18 17 13 15 12 18 15 17 "
        "24 16 18 15 17 15 22 17 16 12 19 18 25 15 20 25 16 23 15 19 22 22 28 "
        "21 19 24 21 22 22 18 21 19 21 23 17 16 26 21 18 22 22 42 1 2 2 4 1 2 "
        "4 2 4 4 5 6 7 2 4 4 3 4 5 3 5 7 6 6 7 7 5 9 8 8 7 10 10 8 12 7 11 9 8 "
        "11 8 9 13 12 10 10 18 8 13 16 13 17 13 16 9 19 18 22 22 15 16 19 14 "
        "13 23 21 16 16 26 19 25 19 24 28 17 27 19 25 23 23 23 26 19 25 24 25 "
        "20 18 17 24 23 24 19 20 29 23 23 24 22 0 0 0 0 1 2 2 2 2 3 2 2 3 4 2 "
        "4 3 2 2 3 4 4 3 2 5 6 5 8 8 7 6 6 5 10 7 10 8 9 10 10 11 7 4 9 8 12 9 "
        "12 10 10 14 10 11 12 10 12 12 13 15 15 11 14 18 15 10 16 18 12 14 17 "
        "18 21 14 14 23 13 27 15 19 18 20 15 17 18 25 21 20 21 17 18 16 19 19 "
        "16 18 15 21 16 19 16 42 1 2 1 3 2 0 3 4 6 5 6 6 6 2 4 4 4 6 6 6 5 6 4 "
        "8 5 9 6 10 8 6 8 12 6 7 10 8 8 10 10 12 7 4 10 10 10 9 11 9 12 13 11 "
        "15 10 15 8 14 12 14 18 10 14 14 15 15 19 16 16 13 18 18 22 17 16 20 "
        "12 20 12 19 16 24 20 14 15 24 14 19 19 20 20 22 22 17 21 14 24 22 19 "
        "19 17 0 1 1 0 3 3 1 1 4 2 2 3 4 2 2 3 4 2 4 3 0 2 3 2 5 6 4 4 8 8 5 9 "
        "5 9 7 11 7 8 11 12 14 6 6 11 8 13 15 13 12 8 11 13 12 14 12 8 14 9 18 "
        "14 16 15 16 15 10 19 15 13 13 18 14 18 11 19 15 10 21 14 16 23 21 19 "
        "19 18 20 20 19 20 18 15 18 21 20 16 19 21 19 15 23 19 0 1 0 0 1 2 0 1 "
        "2 3 3 2 3 2 2 3 5 3 5 4 3 6 6 4 8 8 10 9 11 5 9 8 7 11 6 13 10 13 14 "
        "15 13 9 8 13 9 13 10 14 10 14 16 15 19 11 16 10 11 18 20 22 16 17 21 "
        "17 13 22 18 20 16 24 21 21 22 20 24 19 27 16 22 22 26 24 20 22 32 20 "
        "26 29 23 24 26 27 27 16 26 26 28 23 24 26 0 -1 0 -1 -2 0 -1 -3 1 -1 0 "
        "0 2 -2 2 2 1 1 4 -3 -5 -4 -6 -4 -3 -2 -4 -4 -3 3 -1 2 -5 3 3 -4 0 -2 "
        "-1 1 2 -2 0 -7 1 0 6 -1 -3 0 -2 -3 -3 -1 0 -2 1 -3 0 -4 2 4 -4 0 3 4 "
        "8 4 8 5 7 -4 10 -2 4 2 6 7 2 4 10 2 -1 3 5 1 10 2 10 2 12 2 7 0 8 8 5 "
        "2 7 4 42 0 2 1 3 1 1 1 2 1 3 4 5 6 3 6 4 2 5 3 2 6 5 3 7 7 6 7 7 8 6 "
        "10 8 10 7 10 10 10 14 12 15 7 11 15 9 15 17 9 14 13 14 16 20 16 18 11 "
        "17 14 19 24 20 19 17 18 16 23 16 19 16 23 16 23 14 25 24 12 23 18 24 "
        "19 23 23 25 19 29 23 23 21 22 22 21 27 27 17 20 28 28 17 21 22 0 1 0 "
        "1 3 3 2 2 4 4 4 5 5 4 3 3 3 5 4 5 4 4 6 6 5 8 8 6 10 5 5 5 9 8 4 12 3 "
        "8 8 9 9 8 6 12 9 9 9 13 9 10 13 11 14 13 15 9 14 14 21 21 17 17 19 15 "
        "14 22 19 19 14 21 22 22 15 20 20 14 22 17 20 19 22 20 20 24 24 21 25 "
        "22 23 19 25 21 22 16 23 24 20 23 20 22 0 0 1 1 2 1 1 2 2 3 1 5 3 5 1 "
        "2 1 2 2 4 4 4 4 4 6 5 6 5 8 5 5 8 7 5 7 8 7 9 8 10 10 8 5 10 7 10 9 "
        "11 10 6 11 8 11 12 13 6 14 9 11 13 10 11 9 11 10 13 11 12 16 16 16 13 "
        "8 13 15 11 19 12 14 17 15 16 12 17 20 14 22 20 16 17 17 16 17 14 16 "
        "18 19 17 20 18 0 -1 1 -1 -2 -3 -1 -1 -2 -1 0 -2 0 2 0 2 1 -1 0 -2 1 0 "
        "-1 -3 -1 -4 -2 -1 -3 1 -1 1 0 1 4 -5 7 2 2 0 2 -1 1 -1 0 1 0 -2 1 3 1 "
        "2 -2 -3 -3 -1 -1 -2 -9 -2 -5 -6 -7 -2 -2 -8 -6 -8 -2 -1 -7 -5 -1 -2 "
        "-3 1 1 -1 -3 -7 1 -1 -3 -10 2 -2 -7 -2 1 1 -1 4 3 0 -6 -4 -1 -4 -1 -6 "
        "42 0 2 2 3 1 2 4 3 6 3 6 5 8 2 4 4 3 5 3 5 4 6 4 5 6 7 4 9 9 9 7 9 9 "
        "8 12 7 10 8 8 12 8 8 11 11 13 12 14 10 13 15 10 15 14 17 10 16 14 20 "
        "20 15 19 15 16 15 17 23 17 14 19 17 25 15 21 23 17 24 14 28 20 21 21 "
        "23 20 26 22 23 20 25 20 25 25 21 22 20 32 24 19 23 19 42 0 2 2 4 1 2 "
        "2 3 3 3 6 4 6 4 5 1 4 4 4 5 4 6 5 6 5 6 7 9 8 8 7 10 8 6 12 7 7 9 12 "
        "11 8 7 11 9 9 11 11 10 13 10 16 14 11 13 11 15 11 19 15 13 15 16 12 "
        "11 18 15 15 11 16 11 20 12 13 16 13 21 12 17 21 22 22 20 13 25 18 19 "
        "21 20 22 23 18 21 18 18 24 20 21 22 22 42 0 1 2 3 0 2 2 1 1 3 4 4 6 3 "
        "4 0 2 2 1 2 2 4 5 1 4 3 4 3 1 -1 3 7 3 -1 3 0 0 4 1 5 1 5 5 4 3 5 -2 "
        "3 4 3 3 7 3 2 3 8 0 2 8 2 4 7 4 8 7 1 8 2 4 7 9 6 1 11 8 10 9 7 6 10 "
        "11 11 12 7 9 6 9 15 12 17 11 10 7 9 11 8 12 9 10 42 0 1 2 2 1 1 2 2 4 "
        "3 6 5 6 4 5 3 4 5 4 4 4 5 6 6 6 7 6 8 8 7 7 8 7 5 11 6 5 8 8 12 11 6 "
        "10 10 14 12 12 6 9 14 12 17 8 15 11 15 13 21 18 17 12 12 11 13 13 14 "
        "11 6 17 11 18 13 18 21 16 19 17 17 18 18 20 22 17 25 19 24 19 22 24 "
        "20 23 22 15 19 29 20 21 22 22 42 1 2 2 5 2 3 5 3 5 5 6 7 9 3 6 5 4 6 "
        "6 5 6 9 6 8 9 10 8 12 9 9 10 14 12 10 14 11 14 14 11 13 9 9 14 12 15 "
        "14 17 10 15 19 16 16 11 20 12 18 16 18 24 15 19 20 18 17 21 19 16 16 "
        "23 19 28 19 21 23 13 26 18 22 23 27 25 21 18 30 24 26 25 31 26 32 28 "
        "29 20 22 34 27 25 27 27 42 1 1 1 2 1 1 3 2 4 5 4 6 6 2 6 5 2 5 4 4 5 "
        "6 6 6 7 8 5 7 8 4 8 10 12 5 10 7 8 11 13 12 5 6 11 12 9 10 13 10 14 "
        "15 9 18 13 17 10 17 16 19 19 15 19 19 12 11 15 14 12 15 22 17 19 14 "
        "19 21 17 21 16 18 20 23 14 25 20 27 22 23 17 28 24 32 22 27 24 22 27 "
        "24 22 28 21 0 0 1 0 1 1 1 2 2 4 2 2 3 3 2 2 4 3 5 3 3 5 5 3 7 4 5 7 "
        "10 7 7 7 11 7 8 9 10 9 12 12 11 9 6 13 6 13 12 14 11 8 9 14 15 13 14 "
        "10 15 12 19 19 15 12 16 19 10 17 15 13 12 18 12 21 14 13 13 10 22 14 "
        "15 19 17 20 14 20 21 19 22 23 23 19 21 23 22 20 20 24 20 21 25 20 0 0 "
        "1 0 2 2 1 0 4 2 2 5 5 5 4 6 3 4 5 5 3 4 3 4 6 5 4 4 6 10 7 10 8 10 9 "
        "12 11 10 10 11 17 9 7 13 10 13 14 14 10 13 13 16 17 13 16 11 17 8 18 "
        "16 19 16 17 12 11 19 20 18 14 17 19 22 11 19 18 13 22 14 23 20 19 17 "
        "18 18 29 24 24 18 22 20 22 23 27 21 23 27 25 18 23 22 0 0 1 0 2 1 1 1 "
        "3 3 3 4 4 6 5 6 3 4 4 5 5 4 4 6 10 8 9 7 8 8 5 8 9 10 8 13 11 11 13 "
        "12 17 9 7 13 7 15 13 11 13 16 18 19 19 16 21 14 21 14 22 27 21 16 21 "
        "19 15 21 20 17 18 23 18 28 18 21 22 10 26 20 22 22 26 22 20 20 29 21 "
        "24 21 25 29 26 28 26 17 25 30 28 26 24 29 42 0 1 1 2 1 2 3 2 4 5 5 5 "
        "7 2 5 2 2 4 3 4 5 4 4 6 6 6 6 4 7 4 6 8 8 3 8 6 8 9 7 11 4 6 11 6 12 "
        "11 10 11 10 15 7 17 17 16 12 15 14 17 19 16 18 15 16 11 17 18 14 12 "
        "18 13 21 13 18 21 14 20 15 20 21 21 19 22 19 25 19 19 22 21 25 25 25 "
        "20 20 23 30 25 24 25 21 42 0 1 2 2 1 2 3 1 4 2 4 3 5 2 5 5 3 5 4 5 6 "
        "8 5 7 7 8 7 10 9 9 8 9 9 7 17 8 10 12 10 15 10 10 16 14 10 13 16 9 14 "
        "14 18 16 11 16 8 18 15 22 22 18 20 17 16 17 20 20 15 12 16 17 23 14 "
        "19 21 15 23 16 24 21 20 21 22 21 26 24 24 22 26 23 25 29 25 22 21 32 "
        "29 20 26 30 42 -1 1 2 0 -2 1 2 -1 3 0 0 -1 0 0 -1 -1 1 3 -2 1 0 0 1 "
        "-1 -3 -1 2 1 2 6 0 4 0 3 1 4 4 5 4 3 3 2 3 -2 -3 4 9 2 5 4 7 6 5 9 6 "
        "11 8 11 9 9 9 8 10 2 11 9 9 10 14 9 12 9 8 10 8 18 6 12 8 14 10 11 10 "
        "15 15 11 13 14 9 12 5 23 14 14 18 14 17 14 15 0 0 0 0 1 2 2 2 3 5 5 5 "
        "5 6 3 5 3 5 5 6 7 5 5 6 9 7 10 10 9 9 7 10 13 10 7 14 11 9 10 10 16 "
        "10 6 15 8 16 8 16 8 13 17 13 19 13 16 10 14 13 15 19 16 15 16 15 13 "
        "18 19 20 13 22 19 26 15 18 21 15 24 15 23 16 22 22 21 22 30 21 25 26 "
        "27 23 28 26 28 23 24 29 28 20 27 26 42 0 1 1 1 1 1 2 1 4 3 3 4 5 3 5 "
        "5 4 5 6 5 6 5 5 6 5 5 8 9 7 8 5 9 6 7 9 8 8 10 12 13 8 4 12 8 10 10 "
        "11 11 13 12 10 14 10 12 11 14 10 17 17 10 11 16 17 13 13 14 11 13 15 "
        "18 19 13 10 15 12 19 17 19 20 19 15 15 16 25 20 19 16 20 15 22 23 20 "
        "22 15 23 24 23 21 24"};
    std::basic_stringstream<char> stream{input_str};
    std::vector<double> inpt; /*(std::istream_iterator<double>(stream),
                             std::istream_iterator<double>());*/
    double num;
    for (size_t i = 0; i < rank * rank; ++i)
    {
        stream >> num;
        inpt.push_back(num);
    }
    const matrix_t<double> matr(inpt.begin(), inpt.end());

    EXPECT_TRUE(DblCmp::are_eq(matr.calculate_det(), static_cast<double>(42),
                               DblCmp::tolerance<float>));
}

TEST(Matrices, Iterators)
{
    matrix_t<float> m{2, 3, 5, 1, 2, 3, 8, 8, 8};
    for (matrix_t<float>::iterator start = m.begin(), fin = m.end();
         start != fin; ++start)
    {
        *start += 0.5;
    }

    for (matrix_t<float>::iterator start = m.begin(), fin = m.end();
         start != fin; start++)
    {
        *start += 0.5;
    }

    auto start = m.begin();
    EXPECT_TRUE(DblCmp::are_eq(start[0], 3.f));
    EXPECT_TRUE(DblCmp::are_eq(start[1], 4.f));
    EXPECT_TRUE(DblCmp::are_eq(start[2], 6.f));
    EXPECT_TRUE(DblCmp::are_eq(start[3], 2.f));
    EXPECT_TRUE(DblCmp::are_eq(start[4], 3.f));
    EXPECT_TRUE(DblCmp::are_eq(start[5], 4.f));
    EXPECT_TRUE(DblCmp::are_eq(start[6], 9.f));
    EXPECT_TRUE(DblCmp::are_eq(start[7], 9.f));
    EXPECT_TRUE(DblCmp::are_eq(start[8], 9.f));

    auto it1 = m.begin();
    auto it2 = it1 + 5;
    EXPECT_EQ(it2 - it1, 5ul);
    EXPECT_TRUE(it2 > it1);

    static_assert(
        std::is_same_v<std::iterator_traits<decltype(start)>::iterator_category,
                       std::random_access_iterator_tag>);
    static_assert(
        std::is_same_v<std::iterator_traits<decltype(it1)>::iterator_category,
                       std::random_access_iterator_tag>);
    static_assert(
        std::is_same_v<std::iterator_traits<decltype(it2)>::iterator_category,
                       std::random_access_iterator_tag>);
    static_assert(std::random_access_iterator<decltype(m.begin())>);
    static_assert(std::random_access_iterator<decltype(m.end())>);
    static_assert(std::random_access_iterator<decltype(m.cbegin())>);
    static_assert(std::random_access_iterator<decltype(m.cend())>);
}

TEST(Matrices, Operations)
{
    matrix_t<float> m{2, 3, 5, 1, 2, 3, 8, 8, 8};

    m.mul_row(0, 2);
    EXPECT_TRUE(
        std::equal(m.begin(), m.end(),
                   matrix_t<float>{4, 6, 10, 1, 2, 3, 8, 8, 8}.begin()));
    m.div_row(2, 2);
    EXPECT_TRUE(
        std::equal(m.begin(), m.end(),
                   matrix_t<float>{4, 6, 10, 1, 2, 3, 4, 4, 4}.begin()));
    m.add_row(0, 1);
    EXPECT_TRUE(
        std::equal(m.begin(), m.end(),
                   matrix_t<float>{5, 8, 13, 1, 2, 3, 4, 4, 4}.begin()));
    m.add_row(0, 1, -2);
    EXPECT_TRUE(std::equal(m.begin(), m.end(),
                           matrix_t<float>{3, 4, 7, 1, 2, 3, 4, 4, 4}.begin()));
}