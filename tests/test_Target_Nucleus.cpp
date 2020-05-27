#include "gtest/gtest.h"

#include "Target_Nucleus.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestTargetNucleus, TestvMinimalNucleus)
{
	// ARRANGE
	double ER = keV;
	double mDM = 10*GeV;
	double mNucleus = 16 * mNucleon;
	double tol = 1.0e-7;
	// ACT & ASSERT
	ASSERT_NEAR( vMinimal_Nucleus(ER, mDM, mNucleus), 136.756*km/sec, tol);
}
