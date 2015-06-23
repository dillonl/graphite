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
	auto alignmentsList = std::make_shared< gwiz::AlignmentList >(bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr));
	uint32_t count = 0;

	gwiz::IAlignment::SharedPtr alignmentPtr;
	while (alignmentsList->getNextAlignment(alignmentPtr))
	{
		++count;
	}
	// std::cout << "chrom20 count: " << count << std::endl;
}

TEST(BamAlignmentReaderTests, TestLoadAlignmentRegionWithPositions)
{
	std::string path = TEST_BAM_FILE;
	auto bamAlignmentReaderPtr = std::make_shared< gwiz::BamAlignmentReader >(path);
	std::string regionString = "20:1-100";
	auto regionPtr = std::make_shared< gwiz::Region >(regionString);
	auto alignmentsList = std::make_shared< gwiz::AlignmentList >(bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr));
	uint32_t count = 0;

	gwiz::IAlignment::SharedPtr alignmentPtr;
	while (alignmentsList->getNextAlignment(alignmentPtr))
	{
		++count;
	}
	// std::cout << "chrom20 count: " << count << std::endl;
}

#endif //GWIZ_BAMALIGNMENTREADERTESTS_HPP
