#ifndef GWIZ_BAMALIGNMENTREADERTESTS_HPP
#define GWIZ_BAMALIGNMENTREADERTESTS_HPP

#include "core/alignment/AlignmentList.h"
#include "core/alignment/BamAlignmentReader.h"
#include "core/region/Region.h"
#include "config/TestConfig.h"

TEST(BamAlignmentReaderTests, TestLoadAlignmentRegion)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "20";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);

	ASSERT_EQ(alignmentsList.size(), 5940);
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentsRegionWithoutAlignments)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "20:1-100";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 0);
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentAtStartRegionWithPositions)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "20:1-10000000";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 1988);
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentAtStartRegionWithPositionsSmallRegion)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "1:10316100-10318900";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 0);
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentAtMiddleRegionWithPositions)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "1:10000000-20000000";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 1979);
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentAtEndRegionWithPositions)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "1:100000000-600000000";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 7499);
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentOneHundredThousand)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "1:12308541-12408541";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 1977);
}

/*
TEST(BamAlignmentReaderTests, TestLoadAlignmentTMP)
{
	std::string path = "~/data/NA12878.section.bam";
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "1:4000000-4100000";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 49246);
}
*/

TEST(BamAlignmentReaderTests, TestLoadAlignmentGetsAlignmentsInRegions)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "1:100000000-600000000";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	for (auto alignmentPtr : alignmentsList)
	{
		ASSERT_LE(regionPtr->getStartPosition(), alignmentPtr->getPosition() - alignmentPtr->getLength());
		ASSERT_GE(regionPtr->getEndPosition(), alignmentPtr->getPosition());
	}
}


#endif //GWIZ_BAMALIGNMENTREADERTESTS_HPP
