#include "gtest/gtest.h"

#include "DM_Distribution.hpp"

// Headers from libphysica
#include "Natural_Units.hpp"

using namespace obscura;
using namespace libphysica::natural_units;

TEST(TestDMDistribution, TestSHMDefaultConstructor)
{
	// ARRANGE
	Standard_Halo_Model shm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(In_Units(shm.DM_density, GeV / cm / cm / cm), 0.4);
}
