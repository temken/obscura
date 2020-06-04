#include "gtest/gtest.h"

#include "Direct_Detection_Ionization.hpp"

using namespace obscura;

TEST(TestDirectDetectionIonization, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Ionization detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Ionization experiment");
}
