#define private public
#include "matrix.h"
#undef private

#include "gtest/gtest.h"
#include <initializer_list>

using Matrices::const_matrix_t;

TEST(Matrices, Determinant)
{
    const_matrix_t<float> m1{2}, m2{0}, m3{2, 13, -6, 4}, m4{0, 13, -6, 4},
        m5{2, 0, -6, 4}, m6{2, 13, 0, 4}, m7{2, 13, -6, 0}, m8{0, 0, -6, 4},
        m9{0, 13, 0, 4}, m10{0, 13, -6, 0}, m11{2, 0, 0, 4}, m12{2, 0, -6, 0},
        m13{2, 13, 0, 0}, m14{0, 0, 0, 4}, m15{0, 0, -6, 0}, m16{0, 13, 0, 0},
        m17{2, 0, 0, 0}, m18(2), m19{25, 0, -6, -8};

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
}