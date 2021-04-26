#include "gtest/gtest.h"

#include "obscura/Direct_Detection_Migdal.hpp"

using namespace obscura;

TEST(TestDirectDetectionMigdal, TestDefaultConstructor)
{
	// ARRANGE
	DM_Detector_Migdal detector;
	// ACT & ASSERT
	ASSERT_EQ(detector.name, "Migdal scattering experiment");
}
