#include "gtest/gtest.h"

#include <cmath>

#include "obscura/Target_Nucleus.hpp"

#include "libphysica/Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

//1. Kinematic functions
TEST(TestTargetNucleus, TestvMinimalNucleus)
{
	// ARRANGE
	double ER		= keV;
	double mDM		= 10 * GeV;
	double mNucleus = 16 * mNucleon;
	double tol		= 1.0e-7;
	// ACT & ASSERT
	ASSERT_NEAR(vMinimal_Nucleus(ER, mDM, mNucleus), 136.756 * km / sec, tol);
}

TEST(TestTargetNucleus, TestMaximumNuclearRecoilEnergy)
{
	// ARRANGE
	double vDM		= 1.0e-3;
	double mDM		= 10 * GeV;
	double mNucleus = 16 * mNucleon;
	double ER_max	= 2.0 * libphysica::Reduced_Mass(mDM, mNucleus) * libphysica::Reduced_Mass(mDM, mNucleus) * vDM * vDM / mNucleus;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Maximum_Nuclear_Recoil_Energy(vDM, mDM, mNucleus), ER_max);
}

//2. Class for nuclear isotopes.
TEST(TestTargetNucleus, TestIsotopeConstructor)
{
	// ARRANGE
	unsigned int Z = 2;
	unsigned int A = 4;
	// ACT & ASSERT
	ASSERT_EQ(Isotope(Z, A).Z, Z);
	ASSERT_EQ(Isotope(Z, A).A, A);
	ASSERT_DOUBLE_EQ(Isotope(Z, A).abundance, 1.0);
	ASSERT_DOUBLE_EQ(Isotope(Z, A).spin, 0.0);
	ASSERT_DOUBLE_EQ(Isotope(Z, A).sp, 0.0);
	ASSERT_DOUBLE_EQ(Isotope(Z, A).sn, 0.0);
	ASSERT_EQ(Isotope(Z, A).name, "He-4");
}

TEST(TestTargetNucleus, TestThomasFermiRadius)
{
	// ARRANGE
	double tol = 0.01 * Bohr_Radius;
	// ACT & ASSERT
	for(unsigned Z = 1; Z < 50; Z++)
	{
		Isotope iso = Get_Nucleus(Z)[0];
		ASSERT_NEAR(iso.Thomas_Fermi_Radius(), 0.89 / pow(Z, 1.0 / 3.0) * Bohr_Radius, tol);
	}
}

TEST(TestTargetNucleus, TestHelmFormFactor)
{
	// ARRANGE
	Isotope helium(2, 4);
	Isotope xenon(54, 131);
	double q = 100.0 * MeV;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(helium.Helm_Form_Factor(0.0), 1.0);
	ASSERT_LT(helium.Helm_Form_Factor(q), 1.0);
	ASSERT_NEAR(xenon.Helm_Form_Factor(q), 0.322894, 1.0e-4);
}

TEST(TestTargetNucleus, TestPrintSummaryIsotope)
{
	// ARRANGE
	unsigned int Z = 2;
	unsigned int A = 4;
	Isotope helium(Z, A);
	// ACT & ASSERT
	helium.Print_Summary(0);
}

//3. Class for nuclei containing all isotopes occuring in nature
TEST(TestTargetNucleus, TestNucleusConstructor)
{
	// ARRANGE
	Isotope iso(1, 1);
	// ACT & ASSERT
	ASSERT_EQ(Nucleus(iso).name, "H");
	ASSERT_EQ(Nucleus(iso).Number_of_Isotopes(), 1);
}

TEST(TestTargetNucleus, TestNucleusNumberOfIsotopes)
{
	auto nuclei = Import_Nuclear_Data();
	// ACT & ASSERT
	for(auto& element : nuclei)
		ASSERT_EQ(element.Number_of_Isotopes(), element.isotopes.size());
}

TEST(TestTargetNucleus, TestNucleusGetIsotope)
{
	// ARRANGE
	Nucleus oxygen = Get_Nucleus(8);
	// ACT & ASSERT
	ASSERT_EQ(oxygen.Get_Isotope(17).Z, 8);
	ASSERT_EQ(oxygen.Get_Isotope(17).A, 17);
}

TEST(TestTargetNucleus, TestNucleusBrackets)
{
	// ARRANGE
	Nucleus oxygen = Get_Nucleus(8);
	// ACT & ASSERT
	ASSERT_EQ(oxygen[0].A, 16);
	ASSERT_EQ(oxygen[1].A, 17);
	ASSERT_EQ(oxygen[2].A, 18);
}

TEST(TestTargetNucleus, TestNucleusAverageNuclearMass)
{
	// ARRANGE
	Nucleus oxygen = Get_Nucleus(8);
	// ACT & ASSERT
	oxygen.Print_Summary();
	ASSERT_NEAR(oxygen.Average_Nuclear_Mass() / mNucleon, 16.0044, 0.001);
}

TEST(TestTargetNucleus, TestPrintSummaryNucleus)
{
	// ARRANGE
	Nucleus oxygen = Get_Nucleus(8);

	// ACT & ASSERT
	oxygen.Print_Summary();
}

//4. Nuclear data
TEST(TestTargetNucleus, TestImportNuclearData)
{
	// ARRANGE
	// ACT
	auto all_nuclei					= Import_Nuclear_Data();
	unsigned int number_of_isotopes = 0;
	for(auto& nucleus : all_nuclei)
		number_of_isotopes += nucleus.Number_of_Isotopes();
	// ASSERT
	ASSERT_EQ(all_nuclei.size(), 92);
	ASSERT_EQ(number_of_isotopes, 295);
}

TEST(TestTargetNucleus, TestGetIsotope)
{
	// ARRANGE
	unsigned int Z = 65;
	unsigned int A = 159;
	Import_Nuclear_Data();
	// ACT
	Isotope iso = Get_Isotope(Z, A);
	// ASSERT
	ASSERT_EQ(iso.name, "Tb-159");
	ASSERT_EQ(iso.Z, Z);
	ASSERT_EQ(iso.A, A);
	ASSERT_DOUBLE_EQ(iso.abundance, 1.0);
	ASSERT_DOUBLE_EQ(iso.spin, 1.5);
}

TEST(TestTargetNucleus, TestGetNucleusByZ)
{
	// ARRANGE
	unsigned int Z = 28;
	// ACT & ASSERT
	ASSERT_EQ(Get_Nucleus(Z).name, "Ni");
	ASSERT_EQ(Get_Nucleus(Z).Number_of_Isotopes(), 5);
}

TEST(TestTargetNucleus, TestGetNucleusByName)
{
	// ARRANGE
	std::string name = "U";
	// ACT & ASSERT
	ASSERT_EQ(Get_Nucleus(name).name, name);
	ASSERT_EQ(Get_Nucleus(name)[0].Z, 92);
	ASSERT_EQ(Get_Nucleus(name).Number_of_Isotopes(), 3);
}