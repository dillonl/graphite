#ifndef GRAPHITE_BAMHEADERTESTS_HPP
#define GRAPHITE_BAMHEADERTESTS_HPP

#include "core/file/BamFileReader.h"
#include "core/alignment/BamAlignmentReader.h"

// TDOD Curently testing functionality under the BamAlignment file. If the functionality stays under BamAlignment then I need to add these tests to the appropriate test file.

TEST(BamHeaderTests, GetBamHeader)
{
    std::cout << "Test is running!" << std::endl;
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    BamFileReader bamFile(bamPath);
    bamFile.open();
    std::cout << bamFile.getSamHeader() << std::endl;
    bamFile.close();
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
