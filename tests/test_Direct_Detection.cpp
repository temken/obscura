#include "gtest/gtest.h"

#include "Direct_Detection.hpp"

using namespace obscura;

TEST(TestDirectDetection, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "base name");
}
