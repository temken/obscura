#include "gtest/gtest.h"

#include "Target_Electron.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

//1. Kinematic functions
	TEST(TestTargetElectron, TestvMinimalElectrons)
	{
		// ARRANGE
		double mDM = 123*MeV;
		double q = 3*keV;
		double dE = 23*eV;
		double tol = 1.0e-8;
		// ACT & ASSERT
		ASSERT_NEAR( vMinimal_Electrons(q, dE, mDM), 0.00767886, tol);
	}
