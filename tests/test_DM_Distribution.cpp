#include "gtest/gtest.h"

#include "obscura/DM_Distribution.hpp"

#include "libphysica/Integration.hpp"
#include "libphysica/Natural_Units.hpp"

#include "obscura/Astronomy.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDMDistribution, TestSHMDefaultConstructor)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(In_Units(shm.DM_density, GeV / cm / cm / cm), 0.4);
}

TEST(TestDMDistribution, TestSHMConstructor)
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

TEST(TestDMDistribution, TestMinimumDMSpeed)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(shm.Minimum_DM_Speed(), 0.0);
}

TEST(TestDMDistribution, TestMaximumDMSpeed)
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

TEST(TestDMDistribution, TestSHMPDFVelocity)
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

TEST(TestDMDistribution, TestSHMPDFSpeed)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.PDF_Speed(-100 * km / sec), 0.0);
	EXPECT_DOUBLE_EQ(shm.PDF_Speed(1000 * km / sec), 0.0);
}

TEST(TestDMDistribution, TestSHMNormalization)
{
	// ARRANGE
	Standard_Halo_Model shm;
	std::function<double(double)> pdf = [&shm](double v) {
		return shm.PDF_Speed(v);
	};
	double tol = 1.0e-6;
	// ACT & ASSERT
	ASSERT_NEAR(libphysica::Integrate(pdf, shm.Minimum_DM_Speed(), shm.Maximum_DM_Speed(), tol), 1.0, tol);
}

TEST(TestDMDistribution, TestSHMCDFSpeed)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shm.CDF_Speed(-1.0), 0.0);
	EXPECT_DOUBLE_EQ(shm.CDF_Speed(shm.Minimum_DM_Speed()), 0.0);
	EXPECT_NEAR(shm.CDF_Speed(shm.Maximum_DM_Speed()), 1.0, 1.0e-6);
	EXPECT_DOUBLE_EQ(shm.CDF_Speed(1.0), 1.0);
}

TEST(TestDMDistribution, TestTotalFlux)
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

// TEST(TestDMDistribution, TestAverageVelocity)
// {

// }

TEST(TestDMDistribution, TestAverageSpeed)
{
	// ARRANGE
	double rhoDM = 0.3 * GeV / cm / cm / cm;
	double v0	 = 220 * km / sec;
	double vobs	 = 250 * km / sec;
	double vesc	 = 544 * km / sec;
	Standard_Halo_Model shm(rhoDM, v0, vobs, vesc);
	double vMin = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_NEAR(shm.Average_Speed(), 0.00113939, 1.0e-8);
	EXPECT_NEAR(shm.Average_Speed(vMin), 0.00141172, 1.0e-8);
}

TEST(TestDMDistribution, TestGalacticRestFrame)
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

TEST(TestDMDistribution, TestEtaFunction)
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

TEST(TestDMDistribution, TestPrintSummary)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	shm.Print_Summary();
}

TEST(TestDMDistribution, TestSHMSetFunctions)
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

TEST(TestDMDistribution, TestSHMObserverVelocity)
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

TEST(TestDMDistribution, TestImportedSpectrum)
{
	// ARRANGE
	Standard_Halo_Model shm;
	std::string file_path = "SHM_Table.txt";
	shm.Export_PDF_Speed(file_path, 1000);
	double v   = 350 * km / sec;
	double tol = 1.0e-3;
	// ACT
	Imported_DM_Distribution imported_distr(shm.DM_density, file_path);
	// ASSERT
	EXPECT_NEAR(shm.PDF_Speed(v), imported_distr.PDF_Speed(v), tol * shm.PDF_Speed(v));
	EXPECT_NEAR(shm.CDF_Speed(v), imported_distr.CDF_Speed(v), tol);
	EXPECT_NEAR(shm.Eta_Function(v), imported_distr.Eta_Function(v), tol * shm.Eta_Function(v));
}
