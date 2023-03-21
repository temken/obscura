#include "obscura/Target_Atom.hpp"
#include "gtest/gtest.h"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

// 1. Kinematic functions
TEST(TestTargetElectron, TestvMinimalElectrons)
{
	// ARRANGE
	double mDM = 123 * MeV;
	double q   = 3 * keV;
	double dE  = 23 * eV;
	double tol = 1.0e-8;
	// ACT & ASSERT
	ASSERT_NEAR(vMinimal_Electrons(q, dE, mDM), 0.00767886, tol);
}

TEST(TestAtomicElectron, TestConstructor)
{
	// ARRANGE
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 500.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);

	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(Xe_5p.k_min, k_min);
	EXPECT_DOUBLE_EQ(Xe_5p.k_max, k_max);
	EXPECT_DOUBLE_EQ(Xe_5p.q_min, q_min);
	EXPECT_DOUBLE_EQ(Xe_5p.q_max, q_max);
	EXPECT_NEAR(Xe_5p.dlogk, log10(Xe_5p.k_Grid[1] / Xe_5p.k_Grid[0]), 1.0e-10);
	EXPECT_NEAR(Xe_5p.dlogq, log10(Xe_5p.q_Grid[1] / Xe_5p.q_Grid[0]), 1.0e-10);
	EXPECT_EQ(Xe_5p.Nk, Xe_5p.k_Grid.size());
	EXPECT_EQ(Xe_5p.Nq, Xe_5p.q_Grid.size());
	EXPECT_EQ(Xe_5p.n, 5);
	EXPECT_EQ(Xe_5p.l, 1);
	EXPECT_EQ(Xe_5p.name, "Xe_5p");
	EXPECT_DOUBLE_EQ(Xe_5p.binding_energy, 12.4433 * eV);
	EXPECT_EQ(Xe_5p.number_of_secondary_electrons, 0);
}

TEST(TestAtomicElectron, TestResponseFunction)
{
	// ARRANGE
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 500.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);
	double q = 2.0 * keV;
	double E = 10.0 * eV;
	// ACT & ASSERT
	EXPECT_FLOAT_EQ(Xe_5p.Atomic_Response_Function(1, q, E), Xe_5p.Ionization_Form_Factor(q, E));
	for(int response = 1; response <= 1; response++)   // NEEDS TO BE CHANGED to 4 IF MORE RESPONSE FUNCTIONS ARE ADDED
		EXPECT_NE(Xe_5p.Atomic_Response_Function(response, q, E), 0.0);
}

TEST(TestAtomicElectron, TestIonizationFormFactor)
{
	// ARRANGE
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 500.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);
	double q = keV;
	double E = 10.0 * eV;
	// ACT & ASSERT
	EXPECT_GT(Xe_5p.Ionization_Form_Factor(q, E), 0.0);
}

TEST(TestAtomicElectron, TestDipoleApproximation)
{
	// ARRANGE
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 500.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);
	double q0 = 0.5 * keV;
	double E  = 20 * eV;
	double F0 = Xe_5p.Ionization_Form_Factor(q0, E);
	double q  = 0.01 * keV;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(Xe_5p.Ionization_Form_Factor(q, E), q * q / q0 / q0 * F0);
	EXPECT_DOUBLE_EQ(Xe_5p.Atomic_Response_Function(1, q, E), q * q / q0 / q0 * F0);
}

TEST(TestAtomicElectron, TestPrintSummary)
{
	// ARRANGE
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 500.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);

	// ACT & ASSERT
	Xe_5p.Print_Summary();
}

TEST(TestAtom, TestConstructor1)
{
	// ARRANGE
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 500.0 * keV;
	double Eb	 = 12.4433 * eV;
	Atomic_Electron Xe_5p("Xe", 5, 1, Eb, k_min, k_max, q_min, q_max, 0);
	double w = 9 * eV;
	Atom atom(54, w, {Xe_5p});
	// ACT & ASSERT
	ASSERT_EQ(atom.nucleus.name, "Xe");
	ASSERT_EQ(atom.electrons.size(), 1);
	ASSERT_EQ(atom.electrons[0].name, "Xe_5p");
	ASSERT_EQ(atom.W, w);
}

TEST(TestAtom, TestConstructor2)
{
	// ARRANGE
	Atom atom("Ar");
	// ACT & ASSERT
	ASSERT_EQ(atom.nucleus.name, "Ar");
	ASSERT_EQ(atom.electrons.size(), 5);
	ASSERT_EQ(atom.electrons[0].name, "Ar_3p");
	ASSERT_EQ(atom.W, 19.6 * eV);
}

TEST(TestAtom, TestLowestBindingEnergy)
{
	// ARRANGE
	Atom argon("Ar");
	Atom xenon("Xe");
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(argon.Lowest_Binding_Energy(), 16.0824 * eV);
	ASSERT_DOUBLE_EQ(xenon.Lowest_Binding_Energy(), 12.4433 * eV);
}

TEST(TestAtom, TestElectron)
{
	// ARRANGE
	Atom argon("Ar");
	Atom xenon("Xe");
	// ACT & ASSERT
	ASSERT_EQ(argon.Electron(2, 1).name, "Ar_2p");
	ASSERT_EQ(xenon.Electron(4, 0).name, "Xe_4s");
}

TEST(TestAtom, TestPrintSummary)
{
	// ARRANGE
	Atom atom("Xe");
	// ACT & ASSERT
	atom.Print_Summary();
}
