#include "gtest/gtest.h"

#include "DM_Particle_Standard.hpp"

using namespace obscura;

TEST(TestDMParticleStandard, TestDefaultConstructor)
{
	// ARRANGE
	DM_Particle_Standard dm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 10.0);
}
