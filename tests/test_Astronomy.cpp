#include "gtest/gtest.h"

#include "Astronomy.hpp"

using namespace obscura;

TEST(TestAstronomy, TestFractionalDayssinceJ2000)
{
	//ARRANGE
	int day = 1;
	int month = 1;
	int year = 2000;
	double hour = 12;
	// ACT & ASSERT
    ASSERT_DOUBLE_EQ( Fractional_Days_since_J2000(day, month, year, hour),  0.0);
}
