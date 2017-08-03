#include "SAMFileWriter.h"

#include <functional>

namespace graphite
{
    SAMFileWriter::SAMFileWriter (const std::string& path) :
        ASCIIFileWriter(path)
    {
    }

    SAMFileWriter::~SAMFileWriter ()
    {
        if (m_in_out_stream.is_open())
        {
            close();
        }
    }

    // Create an open function that will allow the SAMFileWriter to be read and written.
    bool SAMFileWriter::open ()
    {
		if (m_in_out_stream.is_open()) { return false; }
		//m_in_out_stream.open(m_file_path, std::ios::in | std::ios::out);
		m_in_out_stream.open(m_file_path, std::fstream::in | std::fstream::out | std::fstream::app);
		//m_in_out_stream.open(m_file_path, std::fstream::in | std::fstream::out | std::fstream::ate);
		//this->m_out_stream.open(this->m_file_path);
		return true;
    }
    
    void SAMFileWriter::recordSamLines (IAlignment::SharedPtr alignmentPtr, GSSWMapping::SharedPtr referenceMappingPtr, GSSWMapping::SharedPtr gsswMappingPtr, std::string originalGraphPathHeader, std::string gsswGraphPathHeader, std::string newCigarString) 
    {
        graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(alignmentPtr);
        std::string readGroup;
        if (gsswMappingPtr->getAltCount() > 0)
            readGroup = "ALT";
        else
            readGroup = "REF"; 

        std::string originalSamLine;
        originalSamLine = 
            bamAlignmentPtr->getName()                              + '\t' //  1. QNAME
            + std::to_string(bamAlignmentPtr->getAlignmentFlag())   + '\t' //  2. FLAG
            + originalGraphPathHeader                               + '\t' //  3. RNAME
            // Need to find out why I need to + 1 on the offset.
            + std::to_string(referenceMappingPtr->getOffset())      + '\t' //  4. POS New position.
            + std::to_string(bamAlignmentPtr->getOriginalMapQuality()) + '\t'//  5. MAPQ
            + bamAlignmentPtr->getCigarString()                     + '\t' //  6. New CIGAR string.
            + bamAlignmentPtr->getMateReferenceName()               + '\t' //  7. Place holder for actual value.
            + std::to_string(bamAlignmentPtr->getMatePosition() + 1)+ '\t' //  8. PNEXT +1 because BamTools mate position is 0-based.
            + std::to_string(bamAlignmentPtr->getTemplateLength())  + '\t' //  9. TLEN
            + bamAlignmentPtr->getSequence()                        + '\t' // 10. SEQ
            + bamAlignmentPtr->getFastqQualities()                  + '\t' // 11. QUAL
            + "RG:Z:" + readGroup;                                         // 12. Optional read group field.

        std::string updatedSamLine;
        updatedSamLine = 
            bamAlignmentPtr->getName()                              + '\t' //  1. QNAME
            + std::to_string(bamAlignmentPtr->getAlignmentFlag())   + '\t' //  2. FLAG
            + gsswGraphPathHeader                                   + '\t' //  3. RNAME
            // Need to find out why I need to + 1 on the offset.
            + std::to_string(gsswMappingPtr->getOffset())           + '\t' //  4. POS New position.
            + std::to_string(bamAlignmentPtr->getOriginalMapQuality()) + '\t'//  5. MAPQ
            + newCigarString                                        + '\t' //  6. New CIGAR string.
            + bamAlignmentPtr->getMateReferenceName()               + '\t' //  7. Place holder for actual value.
            + std::to_string(bamAlignmentPtr->getMatePosition() + 1)+ '\t' //  8. PNEXT +1 because BamTools mate position is 0-based.
            + std::to_string(bamAlignmentPtr->getTemplateLength())  + '\t' //  9. TLEN
            + bamAlignmentPtr->getSequence()                        + '\t' // 10. SEQ
            + bamAlignmentPtr->getFastqQualities()                  + '\t' // 11. QUAL
            + "RG:Z:" + readGroup;                                         // 12. Optional read group field.
        
        {
            std::lock_guard< std::mutex > lock(this->m_sam_file_writer_mutex);
            m_sam_alignments.push_back(originalSamLine);
            m_sam_alignments.push_back(updatedSamLine);
        }
    }

    bool SAMFileWriter::write (const char* data, size_t dataLength)
    {
		if (!m_in_out_stream.is_open()) { return false; }
		m_in_out_stream.write(data, dataLength);
        m_in_out_stream << std::endl;
        return true;
    }

    bool SAMFileWriter::writeSamHeader (const char* header, size_t headerLength)
    {
		if (!m_in_out_stream.is_open()) { return false; }
		m_in_out_stream.write(header, headerLength);
        return true;
    }

    bool SAMFileWriter::writeStoredAlignments ()
    {
        if (!m_in_out_stream.is_open()) { return false; }

        if (m_sam_alignments.size() == 0)
            return false;

        for (auto &alignment: m_sam_alignments)
        {
           m_in_out_stream.write(alignment.c_str(), alignment.size());
           m_in_out_stream << std::endl;
        }

        return true;
    }

    void SAMFileWriter::clearStoredAlignments()
    {
        m_sam_alignments.clear();
    }

    bool SAMFileWriter::endOfFile ()
    {
        return m_in_out_stream.eof();
    }

    void SAMFileWriter::setInStreamToBeginning ()
    {
        m_in_out_stream.seekg(0);
    }

    void SAMFileWriter::getInLine (std::string& line)
    {
        std::getline(m_in_out_stream, line);
    }

    void SAMFileWriter::printIosState ()
    {
        std::cout << "Is m_in_out_stream open? " << m_in_out_stream.is_open() << std::endl;
        std::cout << "good()=" << m_in_out_stream.good() << std::endl;
        std::cout << "eof()=" << m_in_out_stream.eof() << std::endl;
        std::cout << "fail()=" << m_in_out_stream.fail() << std::endl;
        std::cout << "bad()=" << m_in_out_stream.bad() << std::endl;
    }

}
