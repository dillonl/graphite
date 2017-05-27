#ifndef GRAPHITE_BAMFILEREADER_HPP
#define GRAPHITE_BAMFILEREADER_HPP

#include "api/BamAlignment.h"
#include "api/BamReader.h"
#include "api/SamHeader.h"

#include <memory>
#include <string>

namespace graphite
{
    /**
     * This class is a wrapper for reading a BAM file using the BamTools api.
     */
    class BamFileReader
    {
    public:
        typedef std::shared_ptr< BamFileReader > SharedPtr;
        BamFileReader (const std::string& bamPath);
        ~BamFileReader ();

        void open ();                   // Open file.
        std::string getSamHeader ();    // Return SAM header as a string.
        void close ();                  // Close file.

        // May want to use the setRegion function from BamTools.
        // See what data is obtained from the BamAlignment
    private:
        std::string m_bam_path;
        BamTools::BamReader m_bam_reader;
        bool m_opened;
    };
}

#endif // GRAPHITE_BAMFILEREADER_HPP
