#include "gtest/gtest.h"

// Headers from libphysica
#include "Natural_Units.hpp"

#include "Astronomy.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

//1. Transform between astronomical coordinate systems.
TEST(TestAstronomy, TransformEquattoGal)
{
	// ARRANGE
	libphysica::Vector ex({1, 0, 0});
	libphysica::Vector ey({0, 1, 0});
	libphysica::Vector ez({0, 0, 1});
	double T				   = 0;
	double tol				   = 1.0e-5;
	std::vector<double> ex_Gal = {-0.0548763, 0.494109, -0.867666};
	std::vector<double> ey_Gal = {-0.873436, -0.444831, -0.198076};
	std::vector<double> ez_Gal = {-0.483836, 0.746982, 0.455984};
	// ACT & ASSERT
	for(unsigned int i = 0; i < 3; i++)
	{
		ASSERT_NEAR(Transform_Equat_to_Gal(ex, T)[i], ex_Gal[i], tol);
		ASSERT_NEAR(Transform_Equat_to_Gal(ey, T)[i], ey_Gal[i], tol);
		ASSERT_NEAR(Transform_Equat_to_Gal(ez, T)[i], ez_Gal[i], tol);
	}
}

TEST(TestAstronomy, TransformGeoEcltoGal)
{
	// ARRANGE
	libphysica::Vector v({1, 2, 3});
	double T = 0;
	// ACT & ASSERT
	for(int i = 0; i < 3; i++)
		ASSERT_DOUBLE_EQ(Transform_GeoEcl_to_Gal(v, T)[i], -1.0 * Transform_HelEcl_to_Gal(v, T)[i]);
}

TEST(TestAstronomy, TransformHelEcltoGal)
{
	// ARRANGE
	libphysica::Vector ex({1, 0, 0});
	libphysica::Vector ey({0, 1, 0});
	double T				   = 0;
	double tol				   = 1.0e-5;
	std::vector<double> ex_Gal = {0.054876, -0.494109, 0.867666};
	std::vector<double> ey_Gal = {0.993824, 0.110992, 0.000352};
	// ACT & ASSERT
	for(unsigned int i = 0; i < 3; i++)
	{
		ASSERT_NEAR(Transform_HelEcl_to_Gal(ex, T)[i], ex_Gal[i], tol);
		ASSERT_NEAR(Transform_HelEcl_to_Gal(ey, T)[i], ey_Gal[i], tol);
	}
}

TEST(TestAstronomy, TransformGaltoEquat)
{
	// ARRANGE
	libphysica::Vector v({1, 2, 3});
	double T   = 0;
	double tol = 1.0e-10;
	// ACT & ASSERT
	for(unsigned int i = 0; i < 3; i++)
		ASSERT_NEAR(Transform_Gal_to_Equat(Transform_Equat_to_Gal(v, T), T)[i], v[i], tol);
}
//2. Times
TEST(TestAstronomy, TestFractionalDayssinceJ2000)
{
	//ARRANGE
	int day		= 1;
	int month	= 1;
	int year	= 2000;
	double hour = 12;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Fractional_Days_since_J2000(day, month, year, hour), 0.0);
	ASSERT_DOUBLE_EQ(Fractional_Days_since_J2000(day, month, year + 1, hour), 366.0);
}

TEST(TestAstronomy, LocalApparentSiderealTime)
{
	//ARRANGE
	int day			 = 26;
	int month		 = 5;
	int year		 = 2020;
	double h		 = 15;
	double m		 = 30;
	double s		 = 48;
	double longitude = 0 * deg;
	double n		 = Fractional_Days_since_J2000(day, month, year, h, m, s);
	double tol		 = 1.0;
	// ACT & ASSERT
	ASSERT_NEAR(In_Units(Local_Apparent_Sidereal_Time(n, longitude), sec), 28167, tol);
}
//3. Sun's and earth's velocity in the galactic frame
TEST(TestAstronomy, SunVelocity)
{
	// ARRANGE
	std::vector<double> vSun = {11.1, 220 + 12.2, 7.3};
	// ACT & ASSERT
	for(int i = 0; i < 3; i++)
		ASSERT_DOUBLE_EQ(In_Units(Sun_Velocity()[i], km / sec), vSun[i]);
}

TEST(TestAstronomy, EarthVelocity)
{
	// ARRANGE
	double n				   = 0;
	std::vector<double> vEarth = {0.0000606364, 0.000727468, 0.000110569};
	double tol				   = 1e-9;
	// ACT & ASSERT
	for(int i = 0; i < 3; i++)
		EXPECT_NEAR(Earth_Velocity(n)[i], vEarth[i], tol);
}
