#include "gtest/gtest.h"

#include "obscura/Target_Crystal.hpp"

#include "libphysica/Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestTargetCrystal, TestSilicon)
{
	// ARRANGE
	Crystal crystal("Si");
	double q = 10 * keV;
	double E = 5 * eV;
	// ACT & ASSERT
	ASSERT_EQ(crystal.name, "Si");
	ASSERT_EQ(crystal.energy_gap, 1.11 * eV);
	ASSERT_GT(crystal.Crystal_Form_Factor(q, E), 0.0);
}

TEST(TestTargetCrystal, TestGermanium)
{
	// ARRANGE
	Crystal crystal("Ge");
	double q = 2.0 * aEM * mElectron;
	double E = 30 * eV;
	// ACT & ASSERT
	ASSERT_EQ(crystal.name, "Ge");
	ASSERT_EQ(crystal.energy_gap, 0.67 * eV);
	ASSERT_GT(crystal.Crystal_Form_Factor(q, E), 0.0);
}
