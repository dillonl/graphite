#ifndef GRAPHITE_BAMALIGNMENTMANAGERTESTS_HPP
#define GRAPHITE_BAMALIGNMENTMANAGERTESTS_HPP

#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignment.h"
#include "api/BamAlignment.h"

#include <memory>
/*
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "20:10000000-30000000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "20:10000000-10003000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "20:1-10000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "1:12300000-12400000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "1:12300000-12400000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "20:30000000-40000000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "1:12300000-12309000";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "20:26151952-26151953";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
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
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::string region2String = "20";
	auto regionPtr2 = std::make_shared< graphite::Region >(region2String, graphite::Region::BASED::ONE);
	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	getAlignmentPtrsFromReader(path, alignmentPtrs, regionPtr2);
	auto alignmentListReaderPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);

	graphite::IAlignmentList::SharedPtr alignmentListManagerPtr;
	getAlignmentPtrsFromManager(path, alignmentListManagerPtr, regionPtr1, regionPtr2);

	compareAlignmentLists(alignmentListReaderPtr, alignmentListManagerPtr);
	ASSERT_EQ(alignmentListReaderPtr->getCount(), 5940);
}
*/

/**
 * Test that the BamAlignment object can be used to output the appropriate columns of the original BAM file.
 */
/*
TEST(BamAlignmentManagerTests, TestWritingBamDataToFile)
{

    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    //std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    std::string regionString = "chr20:1-1000000";
	auto regionPtr = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

    graphite::BamAlignmentReader bamAlignmentReader(bamPath);
    bamAlignmentReader.open();
    auto samplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
    graphite::SampleManager::SharedPtr smPtr = std::make_shared< graphite::SampleManager >(samplePtrs);
    auto iAlignmentPtrs = bamAlignmentReader.loadAlignmentsInRegion(regionPtr, smPtr, true);
    std::vector< graphite::BamAlignment::SharedPtr > bamAlignmentPtrs;

    // Print out first line of BAM region.
    // Will likely want to implement a more robust test method than printing out the values and eyeballing them.
    // May also want to change the function names for these parameters in BamAlignment.h
    for (int i = 0; i < 3; ++i)
    {
        // Cast iAlgnmentPtrs vector to bamAlignment ptrs.
        graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(iAlignmentPtrs[i]);
        std::cout
            << bamAlignmentPtr->getName() << "\t"                   //  1. QNAME
            << bamAlignmentPtr->getRefSeqName()
            << bamAlignmentPtr->getAlignmentFlag() << "\t"          //  2. FLAG
            << "RNAME Place_Holder" << "\t"                         // RNAME ERROR
            //<< bamAlignmentPtr->getRefSeqName() << "\t"           //  3. RNAME ERROR May just need to get the chr from the regionPtr.
            << bamAlignmentPtr->getPosition() + 1 << "\t"           //  4. POS +1 because the BamTools position is 0-based.
            << bamAlignmentPtr->getOriginalMapQuality() << "\t"     //  5. MAPQ
            //<< bamAlignmentPtr->getOriginalCigarData() << "\t"    //  6. CIGAR Need to implment in header file first.
            << "CIGAR_STR Place_holder" << "\t"
            //<< bamAlignmentPtr->getMateID() << "\t"               //  7. RNEXT ERROR Not the correct value
            << "RNEXT Place_Holder" << "\t"                         // RNEXT ERROR Not the correct value
            << bamAlignmentPtr->getMatePosition() + 1 << "\t"       //  8. PNEXT +1 becuase BamTools mate position is 0-based.
            << bamAlignmentPtr->getTemplateLength() << "\t"         //  9. TLEN
            << bamAlignmentPtr->getSequence() << "\t"               // 10. SEQ
            << bamAlignmentPtr->getFastqQualities() << "\t"         // 11. QUAL
            << std::endl;
    }

    bamAlignmentReader.close();

}
*/
/*
TEST(BamAlignmentManagerTests, TestWritingBamDataToFile)
{
    std::cout << "Test is running! " << std::endl;
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    std::string regionString = "chr20:1-1000000";
	auto regionPtr1 = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

	std::vector< graphite::IAlignment::SharedPtr > alignmentPtrs;
	auto bamAlignmentReaderPtr = std::make_shared< graphite::BamAlignmentReader >(bamPath);
	//auto samplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
    graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite::SampleManager >(graphite::BamAlignmentManager::GetSamplePtrs(bamPath));
	alignmentPtrs = bamAlignmentReaderPtr->loadAlignmentsInRegion(regionPtr1, sampleManagerPtr, true);

	while (alignmentListPtrs->getNextAlignment(alignmentPtrs))
    {
        std::cout << alignmentListPtrs.getName() << std::endl;
    }

}
*/


#endif //GRAPHITE_BAMALIGNMENTMANAGERTESTS_HPP
