#include "obscura/Direct_Detection_Ionization.hpp"
#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Experiments.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDirectDetectionIonization, TestConstructor1)
{
	std::string label  = "Ionization experiment";
	double expo		   = kg * day;
	std::string target = "Electrons";
	std::string atom   = "Ar";
	// ARRANGE
	DM_Detector_Ionization detector(label, expo, target, atom);
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Ionization experiment");
	ASSERT_EQ(detector.Target_Particles(), "Electrons");
}

TEST(TestDirectDetectionIonization, TestConstructor2)
{
	// ARRANGE
	std::string label			   = "Ionization experiment";
	double expo					   = kg * day;
	std::string target			   = "Electrons";
	std::vector<std::string> atoms = {"Ar", "Xe"};
	std::vector<double> abund	   = {1.0, 1.0};
	// ACT
	DM_Detector_Ionization detector(label, expo, target, atoms, abund);
	// ASSERT
	ASSERT_EQ(detector.name, "Ionization experiment");
	ASSERT_EQ(detector.Target_Particles(), "Electrons");
}

TEST(TestDirectDetectionIonization, TestMinimumDMSpeed)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector = XENON1T_S2_ER();
	DM_Particle_SI dm(0.5);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(detector.Minimum_DM_Speed(dm), sqrt(2.0 * 12.4433 * eV / dm.mass));
}

TEST(TestDirectDetectionIonization, TestMinimumDMMass)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector = XENON10_S2_ER();
	DM_Particle_SI dm(0.5);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(detector.Minimum_DM_Mass(dm, shm), 2.0 * 12.4433 * eV / shm.Maximum_DM_Speed() / shm.Maximum_DM_Speed());
}

TEST(TestDirectDetectionIonization, TestRS2Bins)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector = XENON10_S2_ER();
	DM_Particle_SI dm(0.5);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	// ACT
	auto spectrum = detector.DM_Signals_Binned(dm, shm);
	// ASSERT
	ASSERT_GT(spectrum.size(), 0);
	for(auto& entry : spectrum)
		ASSERT_GE(entry, 0.0);
}

TEST(TestDirectDetectionIonization, TestRneBins)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector = DarkSide50_S2_ER();
	DM_Particle_SI dm(0.5);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	// ACT
	auto spectrum = detector.DM_Signals_Binned(dm, shm);
	// ASSERT
	ASSERT_GT(spectrum.size(), 0);
	for(auto& entry : spectrum)
		ASSERT_GE(entry, 0.0);
}

TEST(TestDirectDetectionIonization, TestDMSignalsTotal)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector1;
	detector1.Use_Electron_Threshold(3);
	DM_Detector_Ionization_ER detector2 = DarkSide50_S2_ER();
	DM_Detector_Ionization_ER detector3 = XENON1T_S2_ER();
	DM_Particle_SI dm(0.5);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	// ACT & ASSERT
	double rate1 = detector1.DM_Signals_Total(dm, shm);
	double rate2 = detector2.DM_Signals_Total(dm, shm);
	double rate3 = detector3.DM_Signals_Total(dm, shm);
	EXPECT_GT(rate1, 0.0);
	EXPECT_GT(rate2, 0.0);
	EXPECT_GT(rate3, 0.0);
	dm.Set_Interaction_Parameter(2.0 * pb, "Electrons");
	EXPECT_DOUBLE_EQ(detector1.DM_Signals_Total(dm, shm), 2.0 * rate1);
	EXPECT_DOUBLE_EQ(detector2.DM_Signals_Total(dm, shm), 2.0 * rate2);
	EXPECT_DOUBLE_EQ(detector3.DM_Signals_Total(dm, shm), 2.0 * rate3);
}

TEST(TestDirectDetectionIonization, TestdRdE)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector = XENON1T_S2_ER();
	DM_Particle_SI dm(0.5);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	double E = 10.0 * eV;
	// ACT
	double dRdE = detector.dRdE(E, dm, shm);
	dm.Set_Interaction_Parameter(2.0 * pb, "Electrons");
	// ASSERT
	EXPECT_GT(dRdE, 0.0);
	EXPECT_DOUBLE_EQ(detector.dRdE(E, dm, shm), 2.0 * dRdE);
}

TEST(TestDirectDetectionIonization, TestEnergyThreshold)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector;
	detector.Use_Energy_Threshold(10.0 * eV, 50.0 * eV);
	DM_Particle_SI dm(0.9);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	// ACT
	double N_1 = detector.DM_Signals_Total(dm, shm);
	detector.Use_Energy_Threshold(20.0 * eV, 50.0 * eV);
	double N_2 = detector.DM_Signals_Total(dm, shm);
	// ASSERT
	EXPECT_GT(N_1, N_2);
}

// TEST(TestDirectDetectionIonization, TestdRdEIonization)
// {
// 	// ARRANGE
// 	DM_Detector_Ionization_ER detector;
// 	detector.Use_Energy_Threshold(10.0 * eV, 50.0 * eV);
// 	DM_Particle_SI dm(0.9);
// 	dm.Set_Interaction_Parameter(pb, "Electrons");
// 	Standard_Halo_Model shm;
// 	Atom atom("Xe");
// 	double E = 10 * eV;
// 	// ACT & ASSERT
// 	EXPECT_GT(detector.dRdE_Ionization(E, dm, shm, atom), 0.0);
// }

TEST(TestDirectDetectionIonization, TestS2Threshold)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector;
	detector.Use_PE_Threshold(30.0, 5., 30, 100);
	DM_Particle_SI dm(0.9);
	dm.Set_Interaction_Parameter(pb, "Electrons");
	Standard_Halo_Model shm;
	Atom atom("Xe");
	double W = atom.W;
	int S2	 = 35;
	// ACT & ASSERT
	EXPECT_GT(detector.R_S2(S2, dm, shm, W, Get_Nucleus(54), atom[0]), 0.0);
	EXPECT_GT(detector.R_S2(S2, dm, shm, atom), 0.0);
	EXPECT_GT(detector.R_S2(S2, dm, shm), 0.0);
}

TEST(TestDirectDetectionIonization, TestPrintSummary)
{
	// ARRANGE
	DM_Detector_Ionization_ER detector;
	detector.Use_PE_Threshold(30.0, 5., 30, 100);
	// ACT & ASSERT
	detector.Print_Summary();
}