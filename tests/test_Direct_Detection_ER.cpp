#include "gtest/gtest.h"

#include "obscura/Direct_Detection_ER.hpp"

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Target_Atom.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDirectDetectionER, TestdRdEe)
{
	// ARRANGE
	DM_Particle_SI DM(100.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 100.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);
	Atom xenon("Xe");
	double mNucleus = 131.0 * mNucleon;
	double Ee_1		= 10 * eV;
	double Ee_2		= 25 * eV;

	// ACT & ASSERT
	ASSERT_GT(dRdEe_Ionization_ER(Ee_1, DM, shm, mNucleus, Xe_5p), dRdEe_Ionization_ER(Ee_2, DM, shm, mNucleus, Xe_5p));
	ASSERT_GT(dRdEe_Ionization_ER(Ee_1, DM, shm, xenon), dRdEe_Ionization_ER(Ee_2, DM, shm, xenon));
}

TEST(TestDirectDetectionER, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Electron recoil experiment");
}

TEST(TestDirectDetectionER, TestConstructor1)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector("label", kg * day, "Xe");
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "label");
}

TEST(TestDirectDetectionER, TestConstructor2)
{
	// ARRANGE
	std::vector<double> mass_fractions = {0.5, 0.5};
	std::vector<std::string> targets   = {"Xe", "Ar"};
	DM_Detector_Ionization_ER detector("label_2", kg * day, targets, mass_fractions);
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "label_2");
}

TEST(TestDirectDetectionER, TestdRdEIonization)
{
	// ARRANGE
	DM_Particle_SI DM(100.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	double q_min = 1.0 * keV;
	double q_max = 1000.0 * keV;
	double k_min = 0.1 * keV;
	double k_max = 100.0 * keV;
	Atomic_Electron Xe_5p("Xe", 5, 1, 12.4433 * eV, k_min, k_max, q_min, q_max, 0);
	Atom xenon("Xe");
	double E = 10 * eV;
	DM_Detector_Ionization_ER detector;
	Nucleus nucleus = Get_Nucleus(54);
	// ACT & ASSERT
	ASSERT_EQ(detector.dRdE_Ionization(E, DM, shm, nucleus, Xe_5p), dRdEe_Ionization_ER(E, DM, shm, nucleus.Average_Nuclear_Mass(), Xe_5p));
}