#include "gtest/gtest.h"

#include "obscura/Configuration.hpp"

#include "libphysica/Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestConfiguration, TestDefaultConstructor)
{
	// ARRANGE
	obscura::Configuration cfg;
	// ACT & ASSERT
	ASSERT_EQ(cfg.ID, "default");
}

TEST(TestConfiguration, TestReadConfig)
{
	// ARRANGE
	std::string filename = "test.cfg";
	// ACT
	Configuration cfg(filename);
	// ASSERT
	EXPECT_EQ(cfg.ID, "test");
	EXPECT_DOUBLE_EQ(cfg.DM->mass, 10.0);
	EXPECT_DOUBLE_EQ(cfg.DM_distr->DM_density, 0.4 * GeV / cm / cm / cm);
	EXPECT_EQ(cfg.DM_detector->name, "Nuclear recoil");
	EXPECT_DOUBLE_EQ(cfg.constraints_mass_min, 10.0);
	EXPECT_DOUBLE_EQ(cfg.constraints_mass_max, 100.0);
	EXPECT_EQ(cfg.constraints_masses, 10);
	EXPECT_DOUBLE_EQ(cfg.constraints_certainty, 0.95);
}