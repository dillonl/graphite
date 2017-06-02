#ifndef GRAPHITE_BAMHEADERTESTS_HPP
#define GRAPHITE_BAMHEADERTESTS_HPP

#include "core/file/BamHeaderReader.h"
#include "core/alignment/BamAlignmentReader.h"

// TDOD Curently testing functionality under the BamAlignment file. If the functionality stays under BamAlignment then I need to add these tests to the appropriate test file.

TEST(BamHeaderTests, GetBamHeader)
{
    std::cout << "Test is running!" << std::endl;

    // Setup variables.
    std::vector< std::string > graphPathHeader;
    graphPathHeader.push_back("chr20:61537:0");
    std::vector< int > graphPathLength;
    graphPathLength.push_back(439);

    // Load bam.
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    BamHeaderReader bamFile(bamPath);
    
    bamFile.open();
    bamFile.createPathHeaderVector(graphPathHeader, graphPathLength);
    bamFile.addPathHeadersToSamHeader();
    std::string samHeader = bamFile.getModifiedSamHeader();
    bamFile.close();

    // Write modified SAM header to file.
    std::ofstream samFile;
    samFile.open("SamHeader.sam", std::ios::trunc);
    samFile << samHeader;
    samFile.close();
}

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

#endif //GRAPHITE_BAMHEADERTESTS_HPP
