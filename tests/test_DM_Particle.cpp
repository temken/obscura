#include "gtest/gtest.h"

#include "DM_Particle.hpp"

using namespace obscura;

TEST(TestDMParticle, TestDefaultConstructor)
{
	// ARRANGE
	DM_Particle dm;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(dm.mass, 10.0);
}
