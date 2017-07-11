#ifndef GRAPHITE_BAMALIGNMENTTESTS_HPP
#define GRAPHITE_BAMALIGNMENTTESTS_HPP

#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignment.h"
#include "api/BamAlignment.h"

#include <memory>

/**
 * This test ensures that the original BAM entries were obtained and not modified.
 */
TEST(BamAlignmentTests, OriginalBamEntry)
{
    std::string bamPath = "/uufs/chpc.utah.edu/common/home/marth-d1/data/project_bam/hgsvc_from_cram/CHS/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.bam";
    std::string regionString = "chr20:60000-60100";
	auto regionPtr = std::make_shared< graphite::Region >(regionString, graphite::Region::BASED::ONE);

    graphite::BamAlignmentReader bamAlignmentReader(bamPath);
    bamAlignmentReader.open();
    auto samplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
    graphite::SampleManager::SharedPtr smPtr = std::make_shared< graphite::SampleManager >(samplePtrs);
    auto iAlignmentPtrs = bamAlignmentReader.loadAlignmentsInRegion(regionPtr, smPtr, true);
    bamAlignmentReader.close();

    // Will likely want to implement a more robust test method than printing out the values and eyeballing them.
    // May also want to change the function names for these parameters in BamAlignment.h
    std::ofstream samFile;
    samFile.open("BamAlignmentsTests.sam", std::ios::app);

    for (int i = 0; i < iAlignmentPtrs.size(); ++i)
    {
        // Cast iAlgnmentPtrs vector to bamAlignment ptrs.
        graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(iAlignmentPtrs[i]);

        samFile
            << bamAlignmentPtr->getName()           << "\t"         //  1. QNAME
            << bamAlignmentPtr->getAlignmentFlag()  << "\t"         //  2. FLAG
            << bamAlignmentPtr->getReferenceName()  << "\t"         //  3. RNAME 
            << bamAlignmentPtr->getPosition() + 1   << "\t"         //  4. POS +1 because the BamTools position is 0-based.
            << bamAlignmentPtr->getOriginalMapQuality() << "\t"     //  5. MAPQ
            << bamAlignmentPtr->getCigarString()    << "\t"         //  6. CIGAR Need to implment in header file first.
            << bamAlignmentPtr->getMateReferenceName() << "\t"      //  7. RNEXT ERROR Not the correct value
            << bamAlignmentPtr->getMatePosition() + 1 << "\t"       //  8. PNEXT +1 becuase BamTools mate position is 0-based.
            << bamAlignmentPtr->getTemplateLength() << "\t"         //  9. TLEN
            << bamAlignmentPtr->getSequence()       << "\t"         // 10. SEQ
            << bamAlignmentPtr->getFastqQualities() << "\t"         // 11. QUAL
            << std::endl;
    }

    samFile << "MODIFIED OUTPUT BEGINS!" << std::endl; 

    /* for loop appends bamAlignments with the following updated properties:
     *   RNAME
     *   POS
     *   CIGAR
     */
    /* not updated yet
    for (int i = 0; i < iAlignmentPtrs.size(); ++i)
    {
        // Cast iAlgnmentPtrs vector to bamAlignment ptrs.
        graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(iAlignmentPtrs[i]);

        samFile
            << bamAlignmentPtr->getName()           << "\t"         //  1. QNAME
            << bamAlignmentPtr->getAlignmentFlag()  << "\t"         //  2. FLAG
            << bamAlignmentPtr->getReferenceName()  << "\t"         //  3. RNAME 
            << bamAlignmentPtr->getPosition() + 1 << "\t"           //  4. POS +1 because the BamTools position is 0-based.
            << bamAlignmentPtr->getOriginalMapQuality() << "\t"     //  5. MAPQ
            << bamAlignmentPtr->getCigarString() << "\t"            //  6. CIGAR Need to implment in header file first.
            << bamAlignmentPtr->getMateReferenceName() << "\t"      //  7. RNEXT ERROR Not the correct value
            << bamAlignmentPtr->getMatePosition() + 1 << "\t"       //  8. PNEXT +1 becuase BamTools mate position is 0-based.
            << bamAlignmentPtr->getTemplateLength() << "\t"         //  9. TLEN
            << bamAlignmentPtr->getSequence() << "\t"               // 10. SEQ
            << bamAlignmentPtr->getFastqQualities() << "\t"         // 11. QUAL
            << std::endl;
    }
    */

    samFile.close();
}

#endif // GRAPHITE_BAMALIGNMENTTESTS_HPP
