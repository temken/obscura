#include "gtest/gtest.h"

#include "obscura/Direct_Detection_Crystal.hpp"

#include "libphysica/Natural_Units.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Target_Crystal.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDirectDetectionCrystal, TestdRdEe)
{
	// ARRANGE
	DM_Particle_SI DM(100.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	double E_DM_max = DM.mass / 2.0 * shm.Maximum_DM_Speed() * shm.Maximum_DM_Speed();
	double Ee_1		= 10.0 * eV;
	double Ee_2		= 20.0 * eV;

	// ACT & ASSERT
	EXPECT_GT(dRdEe_Crystal(Ee_1, DM, shm, target), dRdEe_Crystal(Ee_2, DM, shm, target));
	EXPECT_FLOAT_EQ(dRdEe_Crystal(2.0 * E_DM_max, DM, shm, target), 0.0);
	double dRdE1 = dRdEe_Crystal(Ee_1, DM, shm, target);
	DM.Set_Interaction_Parameter(1e-37 * cm * cm, "Electrons");
	EXPECT_FLOAT_EQ(dRdEe_Crystal(Ee_1, DM, shm, target), 0.1 * dRdE1);
}

TEST(TestDirectDetectionCrystal, TestRQ)
{
	// ARRANGE
	DM_Particle_SI DM(100.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	int Q1 = 2;
	int Q2 = 4;
	// ACT & ASSERT
	EXPECT_GT(R_Q_Crystal(Q1, DM, shm, target), R_Q_Crystal(Q2, DM, shm, target));
	double RQ1 = R_Q_Crystal(Q1, DM, shm, target);
	DM.Set_Interaction_Parameter(1e-37 * cm * cm, "Electrons");
	EXPECT_FLOAT_EQ(R_Q_Crystal(Q1, DM, shm, target), 0.1 * RQ1);
}

TEST(TestDirectDetectionCrystal, TestRtotal)
{
	// ARRANGE
	DM_Particle_SI DM(500.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	int Q_thr  = 3;
	double sum = 0.0;
	for(int Q = Q_thr; Q < 10; Q++)
		sum += R_Q_Crystal(Q, DM, shm, target);
	// ACT & ASSERT
	EXPECT_NEAR(R_total_Crystal(Q_thr, DM, shm, target), sum, 1e-3 * sum);
}

TEST(TestDirectDetectionCrystal, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Crystal detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Crystal experiment");
}

TEST(TestDirectDetectionCrystal, TestConstructor)
{
	// ARRANGE
	DM_Detector_Crystal detector("Label", 100 * gram * day, "Si");
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Label");
}

TEST(TestDirectDetectionCrystal, TestMinimumDMSpeed)
{
	// ARRANGE
	DM_Particle_SI DM(500.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	DM_Detector_Crystal detector("Label", 100 * gram * day, "Si");
	detector.Use_Q_Threshold(1);
	// ACT & ASSERT
	ASSERT_EQ(detector.Minimum_DM_Speed(DM), sqrt(2.0 * target.energy_gap / DM.mass));
}

TEST(TestDirectDetectionCrystal, MinimumDMMass)
{
	// ARRANGE
	DM_Particle_SI DM(500.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	DM_Detector_Crystal detector("Label", 100 * gram * day, "Si");
	detector.Use_Q_Threshold(1);
	// ACT & ASSERT
	ASSERT_EQ(detector.Minimum_DM_Mass(DM, shm), 2.0 * target.energy_gap / shm.Maximum_DM_Speed() / shm.Maximum_DM_Speed());
}

TEST(TestDirectDetectionCrystal, TestdRdE)
{
	// ARRANGE
	DM_Particle_SI DM(500.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	DM_Detector_Crystal detector("Label", 100 * gram * day, "Si");
	double Ee = 10.0 * eV;
	// ACT & ASSERT
	ASSERT_EQ(detector.dRdE(Ee, DM, shm), dRdEe_Crystal(Ee, DM, shm, target));
}

TEST(TestDirectDetectionCrystal, TestDMSignalsTotal)
{
	// ARRANGE
	DM_Particle_SI DM(500.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	DM_Detector_Crystal detector("Label", 100 * gram * day, "Si");
	detector.Use_Q_Threshold(1);

	// ACT & ASSERT
	ASSERT_EQ(detector.DM_Signals_Total(DM, shm), 100.0 * gram * day * R_total_Crystal(1, DM, shm, target));
}

TEST(TestDirectDetectionCrystal, TestDMSignalsBinned)
{
	// ARRANGE
	DM_Particle_SI DM(500.0 * MeV);
	DM.Set_Interaction_Parameter(1e-36 * cm * cm, "Electrons");
	Standard_Halo_Model shm;
	Crystal target("Si");
	DM_Detector_Crystal detector("Label", 100 * gram * day, "Si");
	detector.Use_Q_Bins(1);
	// ACT & ASSERT
	for(int i = 0; i < 10; i++)
		ASSERT_EQ(detector.DM_Signals_Binned(DM, shm)[i], 100 * gram * day * R_Q_Crystal(i + 1, DM, shm, target));
}

TEST(TestDirectDetectionCrystal, TestPrintSummary)
{
	// ARRANGE
	DM_Detector_Crystal detector;
	// ACT & ASSERT
	detector.Print_Summary();
}