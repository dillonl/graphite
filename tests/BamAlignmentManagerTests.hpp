#ifndef GRAPHITE_BAMALIGNMENTMANAGERTESTS_HPP
#define GRAPHITE_BAMALIGNMENTMANAGERTESTS_HPP

#include "core/alignment/BamAlignmentManager.h"

void getAlignmentPtrsFromReader(const std::string& path, std::vector< graphite::IAlignment::SharedPtr >& alignmentPtrs, graphite::Region::SharedPtr regionPtr)
{
	auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto samplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(path);
	for (auto samplePtr : samplePtrs)
	{
		graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
	}
	alignmentPtrs = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
}

void getAlignmentPtrsFromManager(const std::string& path, graphite::IAlignmentList::SharedPtr& alignmentListPtr, graphite::Region::SharedPtr regionPtr1, graphite::Region::SharedPtr regionPtr2)
{
	//Add this back in once the bamalignmentmanager is working again
	auto samplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(path);
	for (auto samplePtr : samplePtrs)
	{
		graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
	}
	std::string vcfPath = TEST_VCF_FILE;
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPath, regionPtr1, nullptr, 3000);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory

	auto bamAlignmentManagerPtr = std::make_shared< graphite::BamAlignmentManager >(samplePtrs, regionPtr1);

	bamAlignmentManagerPtr->asyncLoadAlignments(variantManagerPtr, 3000);
	bamAlignmentManagerPtr->waitForAlignmentsToLoad();
	bamAlignmentManagerPtr->releaseResources();

	alignmentListPtr = bamAlignmentManagerPtr->getAlignmentsInRegion(regionPtr2);
}

void compareAlignmentLists(graphite::IAlignmentList::SharedPtr alignmentListPtr1, graphite::IAlignmentList::SharedPtr alignmentListPtr2)
{
	graphite::IAlignment::SharedPtr alignmentPtr1;
	graphite::IAlignment::SharedPtr alignmentPtr2;
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "20:10000000-30000000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentSmallRegion)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "20:10000000-10003000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentRegionEmpty)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "20:1-10000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
	ASSERT_EQ(alignmentListManagerPtr->getCount(), 0);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentSmallHundredThousandRegion)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "1";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "1:12300000-12400000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
	ASSERT_EQ(alignmentListManagerPtr->getCount(), 1977);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentSmallHundredThousandExactRegion)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "1:12300000-12400000";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "1:12300000-12400000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
	ASSERT_EQ(alignmentListManagerPtr->getCount(), 1977);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentRegionTwoSpecific)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20:10000000-50000000";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "20:30000000-40000000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
}

TEST(BamAlignmentManagerTests, TestLoadNineThousandRegion)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "1:12300000-12309000";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "1:12300000-12309000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
	ASSERT_EQ(alignmentListReaderPtr->getCount(), 13);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentRegionOverlapByOne)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20:26000000-27000000";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "20:26151952-26151953";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
}

TEST(BamAlignmentManagerTests, TestLoadAlignmentAllChrom20)
{
	std::string path = TEST_BAM_FILE;
	std::string regionString = "20";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString);

	std::string region2String = "20";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
	ASSERT_EQ(alignmentListReaderPtr->getCount(), 5940);
}

#endif //GRAPHITE_BAMALIGNMENTMANAGERTESTS_HPP
