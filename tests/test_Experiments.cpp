#include "gtest/gtest.h"

#include "obscura/Experiments.hpp"

using namespace obscura;

//1. Nuclear recoil experiments
TEST(TestExperiments, TestDamicN)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = DAMIC_N_2011();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DAMIC_N_2011");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestXENON1TN)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = XENON1T_N_2017();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON1T_N_2017");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestCRESSTII)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = CRESST_II();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CRESST-II");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestCRESSTIII)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = CRESST_III();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CRESST-III");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestCRESSTsurface)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = CRESST_surface();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CRESST-surface");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

//2. Electron recoil experiments - Ionization
TEST(TestExperiments, TestXENON10S2ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = XENON10_S2_ER();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON10_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestXENON100S2ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = XENON100_S2_ER();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON100_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestXENON1TS2ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = XENON1T_S2_ER();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON1T_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestDS50ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = DarkSide50_S2_ER();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DarkSide-50_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

//3. Electron recoil experiments - Semiconductor/crystals
TEST(TestExperiments, TestProtoSENSEIsurface)
{
	//ARRANGE
	DM_Detector_Crystal experiment = protoSENSEI_at_Surface();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "protoSENSEI@surface");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestProtoSENSEIminos)
{
	//ARRANGE
	DM_Detector_Crystal experiment = protoSENSEI_at_MINOS();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "protoSENSEI@MINOS");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestSENSEIminos)
{
	//ARRANGE
	DM_Detector_Crystal experiment = SENSEI_at_MINOS();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "SENSEI@MINOS");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestCDMS2018)
{
	//ARRANGE
	DM_Detector_Crystal experiment = CDMS_HVeV_2018();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CDMS-HVeV_2018");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

TEST(TestExperiments, TestCDMS2020)
{
	//ARRANGE
	DM_Detector_Crystal experiment = CDMS_HVeV_2020();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CDMS-HVeV_2020");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
}

//4. Migdal experiments - Ionization
TEST(TestExperiments, TestXENON10S2Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = XENON10_S2_Migdal();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON10_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestXENON100S2Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = XENON100_S2_Migdal();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON100_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestXENON1TS2Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = XENON1T_S2_Migdal();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON1T_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}

TEST(TestExperiments, TestDS50Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = DarkSide50_S2_Migdal();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DarkSide-50_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
}