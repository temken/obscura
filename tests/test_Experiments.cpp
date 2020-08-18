#include "gtest/gtest.h"

#include "Experiments.hpp"

using namespace obscura;

TEST(TestExperiments, TestDamicN)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = DAMIC_N();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DAMIC-2012");
}