#include "gtest/gtest.h"

#include "Direct_Detection_Semiconductor.hpp"

using namespace obscura;

TEST(TestDirectDetectionSemiconductor, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Semiconductor detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Semiconductor experiment");
}
