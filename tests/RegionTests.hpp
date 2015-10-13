#ifndef GRAPHITE_TESTS_REGION_HPP
#define GRAPHITE_TESTS_REGION_HPP

#include <stdexcept>

#include "core/region/Region.h"

TEST(RegionTest, ConstructorWithGetRegionString)
{
	std::string regionString = "20:10-100";
	graphite::Region region(regionString);
	EXPECT_STREQ(regionString.c_str(), region.getRegionString().c_str());
	std::string regionNotEqual = "1:1-1";
	EXPECT_STRNE(regionNotEqual.c_str(), region.getRegionString().c_str());
}


TEST(RegionTest, GetReferenceID)
{
	std::string reference = "20";
	std::string regionString = reference + ":10-100";
	graphite::Region region(regionString);
	EXPECT_STREQ(reference.c_str(), region.getReferenceID().c_str());
	std::string referenceNotEqual = "1";
	regionString = referenceNotEqual + ":10-100";
	EXPECT_STRNE(referenceNotEqual.c_str(), region.getReferenceID().c_str());
}

/*
TEST(RegionTest, GetReferenceIDNormalized)
{
	std::string reference = "chr20";
	std::string regionString = reference + ":10-100";
	graphite::Region region(regionString);
	EXPECT_STREQ("20", region.getReferenceIDNormalized().c_str());
	std::string referenceNotEqual = "1";
	regionString = referenceNotEqual + ":10-100";
	EXPECT_STRNE(reference.c_str(), region.getReferenceIDNormalized().c_str());
}
*/

TEST(RegionTest, GetStartPosition)
{
	graphite::position startPosition = 2000;
	std::string startPositionString = std::to_string(startPosition);
	std::string regionString = "20:" + startPositionString + "-100000";
	graphite::Region region(regionString);
	EXPECT_EQ(startPosition, region.getStartPosition());
	graphite::position startPositionNotEqual = 1;
	std::string startPositionNotEqualString = std::to_string(startPositionNotEqual);
	EXPECT_NE(startPositionNotEqual, region.getStartPosition());

}

TEST(RegionTest, TestRegionRegionIDOnly)
{
	bool success = false;
	std::string regionString = "Y";

	graphite::Region region(regionString);
	EXPECT_EQ(regionString, region.getReferenceID());
	graphite::position startPositionNotEqual = 1;
	std::string startPositionNotEqualString = std::to_string(startPositionNotEqual);
	EXPECT_EQ(region.getStartPosition(), 0);
	EXPECT_EQ(region.getEndPosition(), graphite::MAX_POSITION);
}

TEST(RegionTest, TestRegionInvalidArgument)
{
	bool success = false;
	std::string invalidString = "";
	ASSERT_THROW(graphite::Region region(invalidString), std::invalid_argument);
}

// Tests invalid position
TEST(RegionTest, TestRegionInvalidStartAndEndPositions)
{
	bool success = false;
	std::string invalidString = "10:10000-1";
	ASSERT_THROW(graphite::Region region(invalidString), std::invalid_argument);
}

#endif //GRAPHITE_TESTS_REGION_HPP
