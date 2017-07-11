#ifndef GRAPHITE_BAMHEADERREADERTESTS_HPP
#define GRAPHITE_BAMHEADERREADERTESTS_HPP

#include <iostream>

#include "core/file/BamHeaderReader.h"

/**
 * When running this test then you may want to pipe it to less -S.
 */
TEST(BamHeaderReaderTests, AddReadGroupsToSamHeader)
{
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    
    graphite::BamHeaderReader bamReader = graphite::BamHeaderReader(bamPath);
    bamReader.open();
    bamReader.addReadGroupsToSamHeader();
    std::cout << bamReader.getModifiedSamHeader();
    bamReader.close();
}

#endif // GRAPHITE_BAMHEADERREADERTESTS_HPP
