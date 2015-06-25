#ifndef GWIZ_BAMALIGNMENTMANAGERTESTS_HPP
#define GWIZ_BAMALIGNMENTMANAGERTESTS_HPP

#include "core/alignment/BamAlignmentManager.h"

void getAlignmentPtrsFromReader(const std::string& path, std::vector< gwiz::IAlignment::SharedPtr >& alignmentPtrs, gwiz::Region::SharedPtr regionPtr)
{
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	alignmentPtrs = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
}

void getAlignmentPtrsFromManager(const std::string& path, gwiz::IAlignmentList::SharedPtr& alignmentListPtr, gwiz::Region::SharedPtr regionPtr1, gwiz::Region::SharedPtr regionPtr2)
{
	auto bamAlignmentManagerPtr = std::make_shared< gwiz::BamAlignmentManager >(path, regionPtr1);

	bamAlignmentManagerPtr->asyncLoadAlignments();
	bamAlignmentManagerPtr->waitForAlignmentsToLoad();
	bamAlignmentManagerPtr->releaseResources();

	alignmentListPtr = bamAlignmentManagerPtr->getAlignmentsInRegion(regionPtr2);
}

void compareAlignmentLists(gwiz::IAlignmentList::SharedPtr alignmentListPtr1, gwiz::IAlignmentList::SharedPtr alignmentListPtr2)
{
	gwiz::IAlignment::SharedPtr alignmentPtr1;
	gwiz::IAlignment::SharedPtr alignmentPtr2;
	while (alignmentListPtr1->getNextAlignment(alignmentPtr1) && alignmentListPtr2->getNextAlignment(alignmentPtr2))
	{
		ASSERT_STREQ(alignmentPtr1->getID().c_str(), alignmentPtr2->getID().c_str());
	}
	ASSERT_EQ(alignmentListPtr1->getCount(), alignmentListPtr2->getCount());
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentRegion)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20";
	auto regionPtr1 = std::make_shared< gwiz::Region >(regionString);

	std::string region2String = "20:10000000-30000000";
	auto regionPtr2 = std::make_shared< gwiz::Region >(region2String);
	std::vector< gwiz::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< gwiz::AlignmentList >(alignmentPtrs);

	gwiz::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
}

#endif //GWIZ_BAMALIGNMENTMANAGERTESTS_HPP
