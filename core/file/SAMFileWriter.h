#ifndef GRAPHITE_SAMFILEWRITER_H
#define GRAPHITE_SAMFILEWRITER_H

#include "core/alignment/BamAlignment.h"
#include "core/alignment/IAlignment.h"
#include "core/file/ASCIIFileWriter.h"
#include "core/mapping/GSSWMapping.h"

#include <string>
#include <vector>

namespace graphite
{
    class SAMFileWriter : public ASCIIFileWriter
    {
    public:
        typedef std::shared_ptr< SAMFileWriter > SharedPtr;
        SAMFileWriter (const std::string& path);
        ~SAMFileWriter ();

        bool open () override;

        void recordSamLines (IAlignment::SharedPtr alignmentPtr, GSSWMapping::SharedPtr referenceMappingPtr, GSSWMapping::SharedPtr gsswMappingPtr, std::string originalGraphPathHeader, std::string gsswGraphPathHeader, std::string newCigarString);    // Store the original and updated SamLines in m_sam_lines.
        bool write (const char* data, size_t datalength) override;  // Write data to SAMFileWriter. Used to write the SAM header to file.
        bool writeSamLines ();  // Write entries to SAM file.
        void clearSamLines ();  // Remove all SAM lines from  the m_sam_out_lines vector.

        void printIosState ();
        /* OLD
        //void adjustInStreamPosition (long adjustment);
        //void adjustOutStreamPosition (long adjustment);

        //void setInStreamToBeginning ();

        //std::string getInLine ();   // Reads a single line and returns it as a string (not including the end of line character.
        //long getInStreamPosition ();
        //long getOutStreamPosition ();
        */

    private:
        //std::ofstream m_in_out_stream;
        std::mutex m_sam_file_writer_mutex;
        std::vector< std::string >  m_sam_out_lines;
    };
}

#endif
