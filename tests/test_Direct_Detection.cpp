#include "obscura/Direct_Detection.hpp"
#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"
#include "obscura/DM_Particle_Standard.hpp"
#include "obscura/Direct_Detection_Nucleus.hpp"
#include "obscura/Target_Nucleus.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDirectDetection, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "base name");
}

TEST(TestDirectDetection, TestEnergyBins)
{
	// ARRANGE
	auto oxygen = Get_Nucleus(8);
	DM_Particle_SI dm(100.0 * GeV);
	Standard_Halo_Model shm;
	DM_Detector_Nucleus detector("test", kg * year, {oxygen});
	// ACT
	detector.Use_Energy_Bins(2.0 * keV, 10.0 * keV, 5);
	detector.Set_Observed_Events({6, 3, 1, 0, 0});
	detector.Set_Bin_Efficiencies({0.1, 0.5, 1.0, 1.0, 1.0});
	detector.Set_Expected_Background({2, 1, 0, 0, 0});
	detector.Print_Summary();
	double limit_1 = detector.Upper_Limit(dm, shm);
	detector.Set_Observed_Events({12, 6, 2, 1, 1});
	double limit_2 = detector.Upper_Limit(dm, shm);
	// ASSERT
	ASSERT_GT(limit_2, limit_1);
}

TEST(TestDirectDetection, TestPValue)
{
	// ARRANGE
	double CL	= 0.9;
	auto oxygen = Get_Nucleus(8);
	DM_Particle_SI dm(100.0 * GeV);
	Standard_Halo_Model shm;
	DM_Detector_Nucleus detector("test", kg * year, {oxygen});
	detector.Use_Energy_Threshold(1.0 * keV, 20 * keV);
	detector.Print_Summary();
	// ACT
	double limit = detector.Upper_Limit(dm, shm, CL);
	dm.Set_Interaction_Parameter(limit, "Nuclei");
	// ASSERT
	ASSERT_NEAR(detector.P_Value(dm, shm), 1.0 - CL, 1e-5);
	dm.Set_Interaction_Parameter(0.5 * limit, "Nuclei");
	ASSERT_GT(detector.P_Value(dm, shm), 1.0 - CL);
	dm.Set_Interaction_Parameter(2.0 * limit, "Nuclei");
	ASSERT_LT(detector.P_Value(dm, shm), 1.0 - CL);
}

TEST(TestDirectDetection, TestLikelihoods)
{
	// ARRANGE
	auto oxygen = Get_Nucleus(8);
	DM_Particle_SI dm(100.0 * GeV);
	Standard_Halo_Model shm;
	DM_Detector_Nucleus detector("test", kg * year, {oxygen});
	detector.Use_Energy_Threshold(1.0 * keV, 20 * keV);
	detector.Print_Summary();
	// ACT & ASSERT
	ASSERT_GT(detector.Likelihood(dm, shm), 0.0);
	ASSERT_DOUBLE_EQ(detector.Log_Likelihood(dm, shm), log(detector.Likelihood(dm, shm)));
}

TEST(TestDirectDetection, TestLikelihoodScan)
{
	// ARRANGE
	auto oxygen = Get_Nucleus(8);
	DM_Particle_SI dm(100.0 * GeV);
	Standard_Halo_Model shm;
	DM_Detector_Nucleus detector("test", kg * year, {oxygen});
	detector.Use_Energy_Threshold(1.0 * keV, 20 * keV);
	detector.Print_Summary();
	// ACT
	auto masses	   = libphysica::Log_Space(10, 100, 5);
	auto couplings = libphysica::Log_Space(1e-40 * cm * cm, 1e-30 * cm * cm, 10);
	auto grid	   = detector.Log_Likelihood_Scan(dm, shm, masses, couplings);
	// ASSERT
	ASSERT_EQ(grid.size(), 5 * 10);
	int i = 0;
	for(auto& m : masses)
		for(auto& c : couplings)
		{
			dm.Set_Mass(m);
			dm.Set_Interaction_Parameter(c, "Nuclei");
			EXPECT_DOUBLE_EQ(m, grid[i][0]);
			EXPECT_DOUBLE_EQ(c, grid[i][1]);
			EXPECT_DOUBLE_EQ(detector.Log_Likelihood(dm, shm), grid[i++][2]);
		}
}
// auto masses			= libphysica::Log_Space(10.0 * MeV, 1.0, 5);
// auto cross_sections = libphysica::Log_Space(1e-47 * cm * cm, 1e-37 * cm * cm, 10);
// auto llhs			= cfg.DM_detector->Log_Likelihood_Scan(*cfg.DM, *cfg.DM_distr, masses, cross_sections);
// std::cout << llhs.size() << std::endl;
// for(auto& entry : llhs)
// 	std::cout << entry[0] / MeV << "\t" << entry[1] / cm / cm << "\t" << entry[2] << std::endl;