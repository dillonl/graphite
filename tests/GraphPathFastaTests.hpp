#ifndef GRAPHITE_GRAPHPATHREFTESTS_HPP
#define GRAPHITE_GRAPHPATHREFTESTS_HPP

#include "TestConfig.h"
#include "core/alignment/BamAlignmentManager.h"
#include "core/file/SAMFileWriter.h"
#include "core/graph/GraphManager.h"
#include "core/graph/GSSWGraph.h"
#include "core/reference/Reference.h"
#include "core/reference/FastaReference.h"
#include "core/region/Region.h"
#include "gmock/gmock.h"

#include <memory.h>

// Try loading multiple variants in.
// Check the input parameters.
// Try adding more code from graphite.cpp in.

// Will eventually delet GraphpathAlignmentTests.hpp since those tests will be included in this file.

class GraphPathFasta : public ::testing::Test
{
};

// Need to rethink test design.
TEST_F(GraphPathFasta, verifyCorrectHeader)
{
    /*
    // Get relevant files
    std::vector< std::string > vcfPaths = {TEST_VCF_FILE};
    std::vector< std::string > bamPaths = {TEST_BAM_FILE};
    graphite::SAMFileWriter::SharedPtr tempSamFilePtr;
    //std::string regionString = "1:0-2000000";
    std::string regionString = "1";
    //auto readLength = 5;  
    auto readLength = graphite::BamAlignmentManager::GetReadLength(bamPaths);  
    std::cout << "Read length: " << readLength << std::endl;
	graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite::SampleManager >(graphite::BamAlignmentManager::GetSamplePtrs(bamPaths));
    auto alignmentReaderManagerPtr = std::make_shared< graphite::AlignmentReaderManager< graphite::BamAlignmentReader > >(bamPaths, 4);

    // Open files...

    auto regionPtr = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);
    auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

    // Load in variants
    auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr, fastaReferencePtr, readLength);
    variantManagerPtr->asyncLoadVCFs();
    variantManagerPtr->waitForVCFsToLoadAndProcess();
    
    // Load in alignments
    auto bamAlignmentManager = std::make_shared< graphite::BamAlignmentManager >(sampleManagerPtr, regionPtr, alignmentReaderManagerPtr, true);
    bamAlignmentManager->asyncLoadAlignments(variantManagerPtr, 3000);
    bamAlignmentManager->waitForAlignmentsToLoad();

    std::deque< std::shared_ptr< std::future< void > > > variantManagerFutureFunctions;
    for (auto& iter : variantManagerPtr->getVCFReadersAndVariantListsMap())
    {
        auto futureFunct = graphite::ThreadPool::Instance()->enqueue(std::bind(&graphite::IVariantList::processOverlappingAlleles, iter.second));
        variantManagerFutureFunctions.push_back(futureFunct);
    }
    while (!variantManagerFutureFunctions.empty())
    {
        variantManagerFutureFunctions.front()->wait();
        variantManagerFutureFunctions.pop_front();
    }

    auto gsswAdjudicator = std::make_shared< graphite::GSSWAdjudicator >(90, 1, 4, 6, 1);
    auto gsswGraphManager = std::make_shared< graphite::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
    gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), readLength, true, tempSamFilePtr);
    std::unordered_map< std::string, std::string > headerSequenceMap = gsswGraphManager->getHeaderSequenceMap();
    std::vector< std::string > headers;
    for (auto &headerSequence : headerSequenceMap)
    {
        headers.push_back(headerSequence.first);
    }

    std::cout << "Map size: " << headerSequenceMap.size() << std::endl;
    std::cout << "Header size: " << headers.size() << std::endl;
    for (auto& header : headers)
    {
        std::cout << "HEADER" << std::to_string(header) << std::endl;
    }

    EXPECT_THAT(headers, ::testing::Contains("chr1:909434:0_"));
    EXPECT_THAT(headers, ::testing::Contains("chr1:909434:1_"));
    */
}

/*
TEST_F(GraphPathFasta, verifyCorrectSequence)
{
}
*/

#endif // GRAPHITE_GRAPHPATHFASTATESTS_HPP
