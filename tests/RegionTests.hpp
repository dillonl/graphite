#ifndef GWIZ_TESTS_REGION_HPP
#define GWIZ_TESTS_REGION_HPP

#include <stdexcept>

#include "core/Region.h"

// Tests factorial of negative numbers.
TEST(RegionTest, ConstructorWithGetRegionString)
{
	std::string regionString = "20:10-100";
	gwiz::Region region(regionString);
	EXPECT_STREQ(regionString.c_str(), region.getRegionString().c_str());
	std::string regionNotEqual = "1:1-1";
	EXPECT_STRNE(regionNotEqual.c_str(), region.getRegionString().c_str());
}

// Tests factorial of negative numbers.
TEST(RegionTest, GetReferenceID)
{
	std::string reference = "20";
	std::string regionString = reference + ":10-100";
	gwiz::Region region(regionString);
	EXPECT_STREQ(reference.c_str(), region.getReferenceID().c_str());
	std::string referenceNotEqual = "1";
	regionString = referenceNotEqual + ":10-100";
	EXPECT_STRNE(referenceNotEqual.c_str(), region.getReferenceID().c_str());
}

// Tests factorial of negative numbers.
TEST(RegionTest, GetStartPosition)
{
	gwiz::position startPosition = 2000;
	std::string startPositionString = std::to_string(startPosition);
	std::string regionString = "20:" + startPositionString + "-100000";
	gwiz::Region region(regionString);
	EXPECT_EQ(startPosition, region.getStartPosition());
	gwiz::position startPositionNotEqual = 1;
	std::string startPositionNotEqualString = std::to_string(startPositionNotEqual);
	EXPECT_NE(startPositionNotEqual, region.getStartPosition());

}

// Tests factorial of negative numbers.
TEST(RegionTest, TestRegionInvalidArgument)
{
	bool success = false;
	std::string invalidString = "invalid";
	ASSERT_THROW(gwiz::Region region(invalidString), std::invalid_argument);
}

#endif //GWIZ_TESTS_REGION_HPP
