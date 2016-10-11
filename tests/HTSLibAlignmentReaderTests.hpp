#ifndef GRAPHITE_TESTS_HTSLIBALIGNMENTREADER_HPP
#define GRAPHITE_TESTS_HTSLIBALIGNMENTREADER_HPP

#include <stdexcept>

#include "core/alignment/HTSLibAlignmentReader.h"
#include "core/alignment/SampleManager.hpp"
#include "config/TestConfig.h"

// #include "core/alignment/BamAlignmentReader.h"

graphite::HTSLibAlignmentReader::SharedPtr getAlignmentReader(const std::string& path)
{
	auto alignmentReaderPtr = std::make_shared< graphite::HTSLibAlignmentReader >(path);
	auto samplePtrs = graphite::HTSLibAlignmentReader::GetBamReaderSamples(path);
	for (auto samplePtr : samplePtrs)
	{
		graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
	}
	return alignmentReaderPtr;
}

TEST(HTSLibAlignmentReaderTest, GetSamples)
{
	auto samplePtrs = graphite::HTSLibAlignmentReader::GetBamReaderSamples(TEST_BAM_FILE);;
	ASSERT_EQ(samplePtrs.size(), 1);
	ASSERT_STREQ(samplePtrs[0]->getName().c_str(), "NA12878");
	ASSERT_STREQ(samplePtrs[0]->getReadgroup().c_str(), "SRR622461");
	ASSERT_STREQ(samplePtrs[0]->getPath().c_str(), TEST_BAM_FILE);
}

TEST(HTSLibAlignmentReaderTests, TestLoadAlignmentsRegionWithoutAlignments)
{
	auto alignmentReaderPtr = getAlignmentReader(TEST_BAM_FILE);
	alignmentReaderPtr->open();
	std::string regionString = "20:1-100";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = alignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	alignmentReaderPtr->close();
	ASSERT_EQ(alignmentsList.size(), 0);
}

TEST(HTSLibAlignmentReaderTests, TestLoadAlignmentAtStartRegionWithPositions)
{
	// std::string path = TEST_BAM_FILE;
	// auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto alignmentReaderPtr = getAlignmentReader(TEST_BAM_FILE);
	alignmentReaderPtr->open();
	std::string regionString = "20:1-10000000";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = alignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	alignmentReaderPtr->close();
	ASSERT_EQ(alignmentsList.size(), 1988);
}

TEST(HTSLibAlignmentReaderTests, TestLoadAlignmentAtStartRegionWithPositionsSmallRegion)
{
	// std::string path = TEST_BAM_FILE;
	// auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto alignmentReaderPtr = getAlignmentReader(TEST_BAM_FILE);
	alignmentReaderPtr->open();
	std::string regionString = "1:10316100-10318900";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = alignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 0);
}

TEST(HTSLibAlignmentReaderTests, TestLoadAlignmentAtMiddleRegionWithPositions)
{
	// std::string path = TEST_BAM_FILE;
	// auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto alignmentReaderPtr = getAlignmentReader(TEST_BAM_FILE);
	alignmentReaderPtr->open();
	std::string regionString = "1:10000000-20000000";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = alignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 1979);
}

TEST(HTSLibAlignmentReaderTests, TestLoadAlignmentAtEndRegionWithPositions)
{
	// std::string path = TEST_BAM_FILE;
	// auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto alignmentReaderPtr = getAlignmentReader(TEST_BAM_FILE);
	alignmentReaderPtr->open();
	std::string regionString = "1:100000000-600000000";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = alignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 7499);
}

TEST(HTSLibAlignmentReaderTests, TestLoadAlignmentOneHundredThousand)
{
	// std::string path = TEST_BAM_FILE;
	// auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto alignmentReaderPtr = getAlignmentReader(TEST_BAM_FILE);
	alignmentReaderPtr->open();
	std::string regionString = "1:12308541-12408541";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = alignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 1977);
}

/*
TEST(BamAlignmentReaderTests, TestLoadAlignmentTMP)
{
	std::string path = "~/data/NA12878.section.bam";
	auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	std::string regionString = "1:4000000-4100000";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	ASSERT_EQ(alignmentsList.size(), 49246);
}
*/

/*
TEST(BamAlignmentReaderTests, TestLoadAlignmentGetsAlignmentsInRegions)
{
	// std::string path = TEST_BAM_FILE;
	// auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(path);
	auto bamAlignmentReaderPtr = getBamAlignmentReader(TEST_BAM_FILE);
	std::string regionString = "1:100000000-600000000";
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto alignmentsList = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr);
	for (auto alignmentPtr : alignmentsList)
	{
		ASSERT_LE(regionPtr->getStartPosition(), alignmentPtr->getPosition() - alignmentPtr->getLength());
		ASSERT_GE(regionPtr->getEndPosition(), alignmentPtr->getPosition());
	}
}
*/


#endif
