#include "gtest/gtest.h"

#include "obscura/Target_Crystal.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Statistics.hpp"

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

TEST(TestTargetCrystal, TestSiliconIonizationYield)
{
	// ARRANGE
	DM_Particle_SI DM(100 * MeV);
	DM.Set_Sigma_Electron(pb);
	Standard_Halo_Model SHM;
	Crystal silicon("Si");
	std::random_device rd;
	std::mt19937 PRNG(rd());

	double tol = 1e-3;
	// ACT & ASSERT
	for(int i = 0; i < 100; i++)
	{
		double Ee  = libphysica::Sample_Uniform(PRNG, 1.1 * eV, 50 * eV);
		double sum = 0.0;
		for(int Q = 1; Q <= 20; Q++)
			sum += silicon.Ionization_Yield(Ee, Q);
		EXPECT_NEAR(sum, 1.0, tol);
	}
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
