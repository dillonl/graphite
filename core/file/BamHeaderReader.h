#ifndef GRAPHITE_BAMHEADERREADER_HPP
#define GRAPHITE_BAMHEADERREADER_HPP

#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/SamHeader.h"

#include <memory>
#include <string>

namespace graphite
{
    /**
     * This class can read a BAM file and return the associated SAM header. It can also return a modified version of the original SAM header but does not modify the original file.
     */
    class BamHeaderReader
    {
    public:
        typedef std::shared_ptr< BamHeaderReader > SharedPtr;
        BamHeaderReader (const std::string& bamPath);
        ~BamHeaderReader ();

        void open ();                   // Open file.
        void createPathHeaderVector (std::vector< std::string > graphPathHeaders, std::vector< int > graphPathLengths);    // Creates a compatible vector of SamSequences that can be added to the header.      
        void addPathHeadersToSamHeader ();   // Add reference sequences to SAM header.
        void close ();                  // Close file.
        std::string getOriginalSamHeader ();    // Return original SAM header from input BAM file as a string.
        std::string getModifiedSamHeader ();    // Retrun updated SAM header as a string.

    private:
        std::string m_bam_path;
        BamTools::BamReader m_bam_reader;
        bool m_opened;
        BamTools::SamHeader m_sam_header;
        std::vector< BamTools::SamSequence > m_sam_sequences;   // Vector of SamSequence objects containing the graphPathHeader information.
    };
}

#endif // GRAPHITE_BAMHEADERREADER_HPP
