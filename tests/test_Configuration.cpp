#include "gtest/gtest.h"

#include "Configuration.hpp"

using namespace obscura;

TEST(TestConfiguration, TestDefaultConstructor)
{
	//ARRANGE
	obscura::Configuration cfg;
	// ACT & ASSERT
	ASSERT_EQ( cfg.ID, "default" );
}
