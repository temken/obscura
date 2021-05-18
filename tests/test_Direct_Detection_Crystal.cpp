#include "gtest/gtest.h"

#include "obscura/Direct_Detection_Crystal.hpp"

using namespace obscura;

TEST(TestDirectDetectionCrystal, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Crystal detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Crystal experiment");
}
