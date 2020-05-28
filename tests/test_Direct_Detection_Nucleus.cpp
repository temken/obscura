#include "gtest/gtest.h"

#include "Direct_Detection_Nucleus.hpp"

using namespace obscura;

TEST(TestDirectDetectionNucleus, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Nucleus detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Nuclear recoil experiment");
}
