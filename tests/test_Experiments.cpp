#include "obscura/Experiments.hpp"
#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

//1. Nuclear recoil experiments
TEST(TestExperiments, TestDamicN)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = DAMIC_N_2011();
	DM_Particle_SI dm(10.0 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DAMIC_N_2011");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestXENON1TN)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = XENON1T_N_2017();
	DM_Particle_SI dm(10.0 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON1T_N_2017");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestCRESSTII)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = CRESST_II();
	DM_Particle_SI dm(10.0 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CRESST-II");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestCRESSTIII)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = CRESST_III();
	DM_Particle_SI dm(10.0 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CRESST-III");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestCRESSTsurface)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = CRESST_surface();
	DM_Particle_SI dm(10.0 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CRESST-surface");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

//2. Electron recoil experiments - Ionization
TEST(TestExperiments, TestXENON10S2ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = XENON10_S2_ER();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON10_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestXENON100S2ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = XENON100_S2_ER();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON100_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestXENON1TS2ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = XENON1T_S2_ER();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON1T_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestDS50ER)
{
	//ARRANGE
	DM_Detector_Ionization_ER experiment = DarkSide50_S2_ER();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DarkSide-50_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

//3. Electron recoil experiments - Semiconductor/crystals
TEST(TestExperiments, TestProtoSENSEIsurface)
{
	//ARRANGE
	DM_Detector_Crystal experiment = protoSENSEI_at_Surface();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "protoSENSEI@surface");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestProtoSENSEIminos)
{
	//ARRANGE
	DM_Detector_Crystal experiment = protoSENSEI_at_MINOS();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "protoSENSEI@MINOS");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestSENSEIminos)
{
	//ARRANGE
	DM_Detector_Crystal experiment = SENSEI_at_MINOS();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "SENSEI@MINOS");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestCDMS2018)
{
	//ARRANGE
	DM_Detector_Crystal experiment = CDMS_HVeV_2018();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CDMS-HVeV_2018");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestCDMS2020)
{
	//ARRANGE
	DM_Detector_Crystal experiment = CDMS_HVeV_2020();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "CDMS-HVeV_2020");
	ASSERT_EQ(experiment.Target_Particles(), "Electrons");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

//4. Migdal experiments - Ionization
TEST(TestExperiments, TestXENON10S2Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = XENON10_S2_Migdal();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON10_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestXENON100S2Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = XENON100_S2_Migdal();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON100_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestXENON1TS2Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = XENON1T_S2_Migdal();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "XENON1T_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}

TEST(TestExperiments, TestDS50Migdal)
{
	//ARRANGE
	DM_Detector_Ionization_Migdal experiment = DarkSide50_S2_Migdal();
	DM_Particle_SI dm(0.5 * GeV);
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DarkSide-50_S2");
	ASSERT_EQ(experiment.Target_Particles(), "Nuclei");
	ASSERT_GT(In_Units(experiment.Upper_Limit(dm, shm), cm * cm), 1.0e-50);
	ASSERT_GE(experiment.P_Value(dm, shm), 0.0);
}