#include "gtest/gtest.h"

#include "obscura/DM_Distribution.hpp"

#include "libphysica/Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

// TEST(TestDMDistribution, TestImportedSpectrum)
// {
// 	// ARRANGE
// 	Standard_Halo_Model shm;
// 	std::string file_path = "SHM_Table.txt";
// 	shm.Export_PDF_Speed(file_path, 1000);
// 	double v   = 350 * km / sec;
// 	double tol = 1.0e-3;
// 	// ACT
// 	Imported_DM_Distribution imported_distr(shm.DM_density, file_path);
// 	// ASSERT
// 	EXPECT_NEAR(shm.PDF_Speed(v), imported_distr.PDF_Speed(v), tol * shm.PDF_Speed(v));
// 	EXPECT_NEAR(shm.CDF_Speed(v), imported_distr.CDF_Speed(v), tol);
// 	EXPECT_NEAR(shm.Eta_Function(v), imported_distr.Eta_Function(v), tol * shm.Eta_Function(v));
// }
