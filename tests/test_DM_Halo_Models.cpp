#include "gtest/gtest.h"

#include "obscura/DM_Halo_Models.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"

#include "obscura/Astronomy.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

// 1. Standard halo model (SHM)
TEST(TestStandardHaloModel, TestDefaultConstructor)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_TRUE(shm.DD_use_eta_function);
	ASSERT_DOUBLE_EQ(In_Units(shm.DM_density, GeV / cm / cm / cm), 0.4);
	ASSERT_DOUBLE_EQ(In_Units(shm.Maximum_DM_Speed(), km / sec), 776.58);
}

TEST(TestStandardHaloModel, TestConstructor)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 230 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 600 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.DM_density, rhoDM);
}

TEST(TestStandardHaloModel, TestMinimumDMSpeed)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(shm.Minimum_DM_Speed(), 0.0);
}

TEST(TestStandardHaloModel, TestMaximumDMSpeed)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 230 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 600 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.Maximum_DM_Speed(), vobs + vesc);
}

TEST(TestStandardHaloModel, TestPDFVelocity)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	libphysica::Vector vel({0, 100 * km / sec, 0});
	// ACT & ASSERT
	ASSERT_NEAR(shm.PDF_Velocity(vel), 3.64055e7, 1.0e2);
}

TEST(TestStandardHaloModel, TestPDFSpeed)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.PDF_Speed(-100 * km / sec), 0.0);
	EXPECT_DOUBLE_EQ(shm.PDF_Speed(1000 * km / sec), 0.0);
}

TEST(TestStandardHaloModel, TestNormalization)
{
	// ARRANGE
	Standard_Halo_Model shm;
	double tol = 1.0e-6;
	// ACT & ASSERT
	ASSERT_NEAR(shm.PDF_Norm(), 1.0, tol);
}

TEST(TestStandardHaloModel, TestCDFSpeed)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.CDF_Speed(-1.0), 0.0);
	EXPECT_DOUBLE_EQ(shm.CDF_Speed(shm.Minimum_DM_Speed()), 0.0);
	EXPECT_NEAR(shm.CDF_Speed(shm.Maximum_DM_Speed()), 1.0, 1.0e-6);
	EXPECT_DOUBLE_EQ(shm.CDF_Speed(1.0), 1.0);
}

TEST(TestStandardHaloModel, TestTotalFlux)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 232 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	double mDM		  = 1.0 * GeV;
	double Flux_total = 1.0 / mDM * 9.8855e6 / cm / cm / sec;
	double tol		  = 1e-5 * Flux_total;
	// ACT & ASSERT
	ASSERT_NEAR(shm.Total_DM_Flux(mDM), Flux_total, tol);
}

TEST(TestStandardHaloModel, TestAverageVelocity)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 232 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	shm.Set_Observer_Velocity(27, 5, 2021, 12, 56);
	libphysica::Vector vel_obs = shm.Get_Observer_Velocity();
	double tol				   = 1.0e-7;
	// ACT & ASSERT
	for(int i = 0; i < 3; i++)
		ASSERT_NEAR(shm.Average_Velocity()[i], -1.0 * vel_obs[i], tol);
}

TEST(TestStandardHaloModel, TestAverageSpeed)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_NEAR(shm.Average_Speed(), 0.00113939, 1.0e-7);
	EXPECT_NEAR(shm.Average_Speed(vMin), 0.00141172, 1.0e-7);
}

TEST(TestStandardHaloModel, TestGalacticRestFrame)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 0 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	double tol	= 1.0e-3 * km / sec;
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.Minimum_DM_Speed(), 0);
	EXPECT_DOUBLE_EQ(shm.Maximum_DM_Speed(), 544 * km / sec);
	EXPECT_NEAR(shm.CDF_Speed(vMin), 0.71127388958876704, 1.0e-6);
	EXPECT_NEAR(shm.Average_Speed(), 245.97191 * km / sec, tol);
	EXPECT_DOUBLE_EQ(shm.Eta_Function(vMin), 0.00079276427680802335156063599536891092516647194737494989322329020 * sec / km);
}

TEST(TestStandardHaloModel, TestEtaFunction)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.Eta_Function(shm.Minimum_DM_Speed()), 1073.3369611520407);
	EXPECT_DOUBLE_EQ(shm.Eta_Function(vMin), 447.76034419713295);
	EXPECT_DOUBLE_EQ(shm.Eta_Function(shm.Maximum_DM_Speed()), 0.0);
	EXPECT_DOUBLE_EQ(shm.Eta_Function(1.0), 0.0);
}

TEST(TestStandardHaloModel, TestPrintSummary)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	shm.Print_Summary();
}

TEST(TestStandardHaloModel, TestSetFunctions)
{
	// ARRANGE
	Standard_Halo_Model shm;
	double v0	= 200 * km / sec;
	double vesc = 800 * km / sec;
	libphysica::Vector vel_obs({0, 400 * km / sec, 0});
	// ACT
	shm.Set_Speed_Dispersion(v0);
	shm.Set_Escape_Velocity(vesc);
	shm.Set_Observer_Velocity(vel_obs);
	// ASSERT
	EXPECT_DOUBLE_EQ(shm.Maximum_DM_Speed(), vesc + vel_obs.Norm());
	EXPECT_NEAR(shm.Average_Speed(), 0.00150091, 1.0e-8);
}

TEST(TestStandardHaloModel, TestObserverVelocity)
{
	// ARRANGE
	double rhoDM				 = 0.3 * GeV / cm / cm / cm;
	double v0					 = 220 * km / sec;
	double vobs					 = 250 * km / sec;
	double vesc					 = 544 * km / sec;
	int day						 = 15;
	int month					 = 3;
	int year					 = 2020;
	libphysica::Vector vel_Earth = Earth_Velocity(Fractional_Days_since_J2000(day, month, year));
	double v_Earth				 = vel_Earth.Norm();
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	// ACT
	shm.Set_Observer_Velocity(day, month, year);
	// ASSERT
	ASSERT_DOUBLE_EQ(shm.Maximum_DM_Speed(), v_Earth + vesc);
}

// 2. Standard halo model++ (SHM++) as proposed by Evans, O'Hare and McCabe [arXiv:1810.11468]
TEST(TestSHMplusplus, TestDefaultConstructor)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	// ACT & ASSERT
	ASSERT_TRUE(shmpp.DD_use_eta_function);
	ASSERT_DOUBLE_EQ(In_Units(shmpp.DM_density, GeV / cm / cm / cm), 0.55);
	ASSERT_DOUBLE_EQ(In_Units(shmpp.Maximum_DM_Speed(), km / sec), 768.0);
}

TEST(TestSHMplusplus, TestConstructor)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 230 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 600 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shmpp.DM_density, rhoDM);
}

TEST(TestSHMplusplus, TestMinimumDMSpeed)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(shmpp.Minimum_DM_Speed(), 0.0);
}

TEST(TestSHMplusplus, TestMaximumDMSpeed)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 230 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 600 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shmpp.Maximum_DM_Speed(), vobs + vesc);
}

TEST(TestSHMplusplus, TestPDFVelocity)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	libphysica::Vector vel({0, 100 * km / sec, 0});
	// ACT & ASSERT
	ASSERT_NEAR(shmpp.PDF_Velocity(vel), 2.9133879137198411e7, 1.0e2);
}

TEST(TestSHMplusplus, TestPDFSpeed)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	// ACT & ASSERT
	EXPECT_GT(shmpp.PDF_Speed(300 * km / sec), 0.0);
	EXPECT_DOUBLE_EQ(shmpp.PDF_Speed(-100 * km / sec), 0.0);
	EXPECT_DOUBLE_EQ(shmpp.PDF_Speed(1000 * km / sec), 0.0);
}

TEST(TestSHMplusplus, TestNormalization)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	double tol = 1e-3;
	// ACT & ASSERT
	ASSERT_NEAR(shmpp.PDF_Norm(), 1.0, tol);
}

TEST(TestSHMplusplus, TestCDFSpeed)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shmpp.CDF_Speed(-1.0), 0.0);
	EXPECT_DOUBLE_EQ(shmpp.CDF_Speed(shmpp.Minimum_DM_Speed()), 0.0);
	EXPECT_NEAR(shmpp.CDF_Speed(shmpp.Maximum_DM_Speed()), 1.0, 1.0e-3);
	EXPECT_DOUBLE_EQ(shmpp.CDF_Speed(1.0), 1.0);
}

TEST(TestSHMplusplus, TestTotalFlux)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 232 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	double mDM		  = 1.0 * GeV;
	double Flux_total = 2.5273619321346239e-45 / mDM;
	double tol		  = 1e-5 * Flux_total;
	// ACT & ASSERT
	ASSERT_NEAR(shmpp.Total_DM_Flux(mDM), Flux_total, tol);
}

TEST(TestSHMplusplus, TestAverageVelocity)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 232 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	shmpp.Set_Observer_Velocity(27, 5, 2021, 12, 56);
	libphysica::Vector vel_obs = shmpp.Get_Observer_Velocity();
	double tol				   = 1.0e-6;
	// ACT & ASSERT
	for(int i = 0; i < 3; i++)
		ASSERT_NEAR(shmpp.Average_Velocity()[i], -1.0 * vel_obs[i], tol);
}

TEST(TestSHMplusplus, TestAverageSpeed)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_NEAR(shmpp.Average_Speed(), 0.0011386489523198436, 1.0e-7);
	EXPECT_NEAR(shmpp.Average_Speed(vMin), 0.0013953880461272795, 1.0e-7);
}

TEST(TestSHMplusplus, TestGalacticRestFrame)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 0 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	double tol	= 1.0e-3 * km / sec;
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shmpp.Minimum_DM_Speed(), 0);
	EXPECT_DOUBLE_EQ(shmpp.Maximum_DM_Speed(), 544 * km / sec);
	EXPECT_NEAR(shmpp.CDF_Speed(vMin), 0.72152699650337748, 1.0e-6);
	EXPECT_NEAR(shmpp.Average_Speed(), 0.00080395484422832609, tol);
	EXPECT_DOUBLE_EQ(shmpp.Eta_Function(vMin), 227.69351987958018);
}

TEST(TestSHMplusplus, TestEtaFunction)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shmpp.Eta_Function(shmpp.Minimum_DM_Speed()), 1055.2029097888699);
	EXPECT_DOUBLE_EQ(shmpp.Eta_Function(vMin), 455.21126825811115);
	EXPECT_NEAR(shmpp.Eta_Function(shmpp.Maximum_DM_Speed()), 0.0, 1e-10);
	EXPECT_DOUBLE_EQ(shmpp.Eta_Function(1.0), 0.0);
}

TEST(TestSHMplusplus, TestPrintSummary)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	// ACT & ASSERT
	shmpp.Print_Summary();
}

TEST(TestSHMplusplus, TestSetFunctions)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	double v0	= 200 * km / sec;
	double vesc = 800 * km / sec;
	libphysica::Vector vel_obs({0, 400 * km / sec, 0});
	// ACT
	shmpp.Set_Speed_Dispersion(v0);
	shmpp.Set_Escape_Velocity(vesc);
	shmpp.Set_Observer_Velocity(vel_obs);
	// ASSERT
	EXPECT_DOUBLE_EQ(shmpp.Maximum_DM_Speed(), vesc + vel_obs.Norm());
	EXPECT_NEAR(shmpp.Average_Speed(), 0.0015072180115424205, 1.0e-8);
}

TEST(TestSHMplusplus, TestObserverVelocity)
{
	// ARRANGE
	double rhoDM				 = 0.3 * GeV / cm / cm / cm;
	double v0					 = 220 * km / sec;
	double vobs					 = 250 * km / sec;
	double vesc					 = 544 * km / sec;
	int day						 = 15;
	int month					 = 3;
	int year					 = 2020;
	libphysica::Vector vel_Earth = Earth_Velocity(Fractional_Days_since_J2000(day, month, year));
	double v_Earth				 = vel_Earth.Norm();
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	// ACT
	shmpp.Set_Observer_Velocity(day, month, year);
	// ASSERT
	ASSERT_DOUBLE_EQ(shmpp.Maximum_DM_Speed(), v_Earth + vesc);
}

TEST(TestSHMplusplus, TestSetBeta)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	double pdf = shmpp.PDF_Speed(200 * km / sec);
	shmpp.Set_Beta(0.3);
	// ACT & ASSERT
	ASSERT_NEAR(shmpp.PDF_Norm(), 1.0, 1e-4);
	ASSERT_NE(shmpp.PDF_Speed(200 * km / sec), pdf);
}

TEST(TestSHMplusplus, TestSHMEquality)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	SHM_Plus_Plus shmpp(rhoDM, v0, vobs, vesc);
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	shmpp.Set_Eta(0.0);
	libphysica::Vector vel({100 * km / sec, 100 * km / sec, 100 * km / sec});
	double v = vel.Norm();
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(shm.PDF_Velocity(vel), shmpp.PDF_Velocity(vel));
	ASSERT_DOUBLE_EQ(shm.PDF_Speed(v), shmpp.PDF_Speed(v));
	ASSERT_DOUBLE_EQ(shm.CDF_Speed(v), shmpp.CDF_Speed(v));
	ASSERT_DOUBLE_EQ(shm.Eta_Function(v), shmpp.Eta_Function(v));
}