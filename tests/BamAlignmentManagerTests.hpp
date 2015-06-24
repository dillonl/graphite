#ifndef GWIZ_BAMALIGNMENTMANAGERTESTS_HPP
#define GWIZ_BAMALIGNMENTMANAGERTESTS_HPP

#include "core/alignment/BamAlignmentManager.h"

TEST(BamAlignmentManagerTests, TestLoadAlignmentRegion)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto bamAlignmentManagerPtr = std::make_shared< gwiz::BamAlignmentManager >(path, regionPtr);

	bamAlignmentManagerPtr->asyncLoadAlignments();
	bamAlignmentManagerPtr->waitForAlignmentsToLoad();
	bamAlignmentManagerPtr->releaseResources();

	std::string region2String = "20:10000000-30000000";
	auto region2Ptr = std::make_shared< gwiz::Region >(region2String);
	auto alignmentsListPtr = bamAlignmentManagerPtr->getAlignmentsInRegion(region2Ptr);

}

#endif //GWIZ_BAMALIGNMENTMANAGERTESTS_HPP
