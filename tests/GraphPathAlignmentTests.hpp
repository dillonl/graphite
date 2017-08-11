#ifndef GRAPHITE_GRAPHPATHALIGNMENTTESTS_HPP
#define GRAPHITE_GRAPHPATHALIGNMENTTESTS_HPP

#include "core/alignment/BamAlignmentReader.h"
#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignment.h"
#include "core/file/BamHeaderReader.h"
#include "core/sample/SampleManager.h"

#include <memory>


// TODO Curently testing functionality under the BamAlignment file. If the functionality stays under BamAlignment then I need to add these tests to the appropriate test file.
// Use the BamAlignment class

TEST(GraphPathAlignmentTests, GetBamHeader)
{
    // Setup variables.
    std::vector< std::string > graphPathHeaders;
    graphPathHeaders.push_back("chr20:61537:0");
    graphPathHeaders.push_back("chr20:61537:1");
    std::vector< int > graphPathLengths;
    graphPathLengths.push_back(439);
    graphPathLengths.push_back(439);

    // Load bam.
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/u0702603/Marthlab/Projects/GraphiteDataViewer/graphite_acmiller015/tests/data/test.bam";
    BamHeaderReader bamFile(bamPath);
    
    bamFile.open();
    bamFile.addPathHeaderToSamHeader(graphPathHeaders[0], graphPathLengths[0]);
    std::string samHeader = bamFile.getModifiedSamHeader();
    bamFile.close();

    // Write modified SAM header to file.
    std::ofstream samFile;
    samFile.open("/uufs/chpc.utah.edu/common/home/u0702603/Marthlab/Projects/GraphiteDataViewer/graphite_acmiller015/tests/data/TestSamHeader.sam", std::ios::trunc);
    samFile << samHeader;
    samFile.close();
}

/** Need to define the purpose of this test before continuing.
 *
TEST(GraphPathAlignmentTests, UpdateBamAlignments)
{
    // Load BAM
    // Get alignments
    // Update alignments
    // write out alignments.
    
    // Setup initial varaibles.
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/u0702603/Marthlab/Projects/GraphiteDataViewer/graphite_acmiller015/tests/data/test.bam";
    std::string fastaPath = "/uufs/chpc.utah.edu/common/home/u0702603/Marthlab/Projects/GraphiteDataViewer/graphite_acmiller015/tests/data/test.fasta";
    std::vector< std::string > bamPaths;
    bamPaths.push_back(bamPath);

    uint32_t readLength = graphite::BamAlignmentManager::GetReadLength(bamPaths);

    auto fastaRegionPtr = std::make_shared< graphite::Region >("1", 1, 1000, graphite::Region::BASED::ONE);
    auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, fastaRegionPtr);

    graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite:: SampleManager >(graphite::BamAlignmentManager::GetSamplePtrs(bamPaths));


    // Load BAM alignments.
    graphite::BamAlignmentReader bamReader(bamPath);
    bamReader.open();
    std::vector< IAlignment::SharedPtr > alignmentPtrs = bamReader.loadAlignmentsInRegion(fastaRegionPtr, sampleManagerPtr, true);
    auto alignmentListPtr = std::make_shared< graphite::AlignmentList >(alignmentPtrs);
    IAlignment::SharedPtr alignmentPtr;

    // Graph Manager
    // don't have the arguments completely defined.
    auto gsswGraphManager = std::make_shared< graphite::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
    gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), readLength);

    /*
    std::string samPath = "/uufs/chpc.utah.edu/common/home/u0702603/Marthlab/Projects/GraphiteDataViewer/graphite_acmiller015/tests/data/TestSamAlignments.sam";;
    while (alignmentListPtr->getNextAlignment(alignmentPtr))
    {
        auto gsswGraphContainer = gsswGraphPtr->getGraphContainer();
        auto tracebackPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr, gsswGraphContainer
        auto gsswMappingPtr = std::make_shared< GSSWMapping >(tracebackPtr, alignmentPtr);


        graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(alignmentPtr);
        std::ofstream samFile;
        samFile.open(samPath, std::ios::app);
        samFile
            << bamAlignmentPtr->getName()           << "\t" //  1. QNAME
            << bamAlignmentPtr->getAlignmentFlag()  << "\t" //  2. FLAG
            // Need to modify the following line
            << "Header_Place_Holder"                << "\t" //  3. RNAME
            << "Pos_Place_holder"                   << "\t" //  4. POS
            << bamAlignmentPtr->getOriginalMapQuality()   <<  "\t"  //  5. MAPQ
            << "Cigar_String_Place_Holder"          << "\t" //  6. New CIGAR string
            << "*"                                  << "\t" //  7. RNEXT
            << bamAlignmentPtr->getMatePosition() + 1 << "\t" //  8. PNEXT +1 because BamTools mate position is 0-based.
            << bamAlignmentPtr->getTemplateLength() << "\t" //  9. TLEN
            << bamAlignmentPtr->getSequence()       << "\t" // 10. SEQ
            << bamAlignmentPtr->getFastqQualities() << "\t" // 11. QUAL
            << std::endl;
    }
    bamReader.close();
    */
/*
}
*/

/*
TEST(BamHeaderTests, GetBamInfo)
{
    std::cout << "Test is running!" << std::endl;
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    BamAlignmentReader bamFile(bamPath);
    bamFile.open();
    std::cout << bamFile.Name << std::endl;
    //std::cout << bamFile.getSamHeader() << std::endl;
    bamFile.close();
}
*/

#endif // GRAPHITE_GRAPHPATHALIGNMENTTESTS_HPP
