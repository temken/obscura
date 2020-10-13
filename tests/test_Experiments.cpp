#include "gtest/gtest.h"

#include "Experiments.hpp"

using namespace obscura;

TEST(TestExperiments, TestDamicN)
{
	//ARRANGE
	DM_Detector_Nucleus experiment = DAMIC_N_2011();
	// ACT & ASSERT
	ASSERT_EQ(experiment.name, "DAMIC_N_2011");
}