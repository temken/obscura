#include "obscura/DM_Distribution.hpp"
#include "gtest/gtest.h"

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Utilities.hpp"

#include "obscura/DM_Halo_Models.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDMDistribution, TestDefaultConstructor)
{
	// ARRANGE
	DM_Distribution dm_distr;
	// ACT & ASSERT
	EXPECT_EQ(dm_distr.DM_density, 0.0);
	EXPECT_FALSE(dm_distr.DD_use_eta_function);
	EXPECT_DOUBLE_EQ(dm_distr.Minimum_DM_Speed(), 0.0);
	EXPECT_DOUBLE_EQ(dm_distr.Maximum_DM_Speed(), 1.0);
}

TEST(TestDMDistribution, TestPDFSpeed)
{
	// ARRANGE
	DM_Distribution dm_distr;
	double v = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(dm_distr.PDF_Speed(v), 0.0);
}

TEST(TestDMDistribution, TestCDFSpeed)
{
	// ARRANGE
	SHM_Plus_Plus shmpp;
	double v = 300 * km / sec;
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(shmpp.CDF_Speed(0.0), 0.0);
	EXPECT_GT(shmpp.CDF_Speed(v), 0.0);
	EXPECT_LT(shmpp.CDF_Speed(v), 1.0);
	EXPECT_DOUBLE_EQ(shmpp.CDF_Speed(shmpp.Maximum_DM_Speed()), 1.0);
}

TEST(TestDMDistribution, TestPDFNorm)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	EXPECT_NEAR(shm.PDF_Norm(), 1.0, 1e-6);
}

TEST(TestDMDistribution, TestEtaFunction)
{
	// ARRANGE
	Standard_Halo_Model shm;
	DM_Distribution dm_distr;
	double v = 1e-3;
	// ACT & ASSERT
	EXPECT_GT(shm.Eta_Function(v), 0.0);
	EXPECT_DOUBLE_EQ(dm_distr.Eta_Function(v), 0.0);
}

TEST(TestDMDistribution, TestPrintSummary)
{
	// ARRANGE
	DM_Distribution dm_distr;
	// ACT & ASSERT
	dm_distr.Print_Summary();
}

TEST(TestDMDistribution, TestExportPDF)
{
	// ARRANGE
	DM_Distribution dm_distr;
	std::string filepath = "test_pdf.txt";
	// ACT
	dm_distr.Export_PDF_Speed(filepath);
	// ASSERT
	EXPECT_TRUE(libphysica::File_Exists(filepath));
	auto pdf = libphysica::Import_Table(filepath);
	EXPECT_EQ(pdf.size(), 100);
	EXPECT_DOUBLE_EQ(pdf[0][0], 0.0);
	EXPECT_NEAR(pdf.back()[0], In_Units(1.0, km / sec), 1.0);
}

TEST(TestDMDistribution, TestExportEta)
{
	// ARRANGE
	DM_Distribution dm_distr;
	std::string filepath = "test_eta.txt";
	// ACT
	dm_distr.Export_Eta_Function(filepath);
	// ASSERT
	EXPECT_TRUE(libphysica::File_Exists(filepath));
	auto eta = libphysica::Import_Table(filepath);
	EXPECT_EQ(eta.size(), 100);
	EXPECT_DOUBLE_EQ(eta[0][0], 0.0);
	EXPECT_NEAR(eta[99][0], In_Units(1.0, km / sec), 1.0);
}

// 2. Import a tabulated DM distribution from a file (format v[km/sec] :: f(v) [sec/km])
TEST(TestImportedDMDistribution, TestImportedSpectrum)
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

TEST(TestImportedDMDistribution, TestPrintSummary)
{
	// ARRANGE
	Standard_Halo_Model shm;
	std::string file_path = "SHM_Table.txt";
	shm.Export_PDF_Speed(file_path, 1000);
	// ACT
	Imported_DM_Distribution imported_distr(shm.DM_density, file_path);
	// ASSERT
	imported_distr.Print_Summary();
}
