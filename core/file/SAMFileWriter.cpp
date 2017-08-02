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
        if (m_opened)
        {
            close();
        }
    }

    // Create an open function that will allow the SAMFileWriter to be read and written.
    bool SAMFileWriter::open ()
    {
		if (m_opened) { return false; }
		//m_out_stream.open(m_file_path, std::fstream::in | std::fstream::out);
		this->m_out_stream.open(this->m_file_path);
		m_opened = true;
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
            m_sam_out_lines.push_back(originalSamLine);
            m_sam_out_lines.push_back(updatedSamLine);
        }
    }

    bool SAMFileWriter::write (const char* data, size_t dataLength)
    {
		if (!m_opened) { return false; }
		this->m_out_stream.write(data, dataLength);
		m_opened = true;
    }

    bool SAMFileWriter::writeSamLines ()
    {
        if (!m_opened) { return false; }

        if (m_sam_out_lines.size() == 0)
            return false;

        for (auto &samLine : m_sam_out_lines)
        {
           m_out_stream.write(samLine.c_str(), samLine.size());
           m_out_stream << std::endl;
        }

        return true;
    }

    void SAMFileWriter::clearSamLines()
    {
        m_sam_out_lines.clear();
    }

    /* OLD
    void SAMFileWriter::adjustInStreamPosition (long adjustment)
    {
        long streamPos = m_out_stream.tellg();
        m_out_stream.seekg(streamPos + adjustment);
    }

    void SAMFileWriter::adjustOutStreamPosition (long adjustment)
    {
        long streamPos = m_out_stream.tellp();
        m_out_stream.seekp(streamPos + adjustment);
    }

    void SAMFileWriter::setInStreamToBeginning ()
    {
        m_out_stream.seekg(0);
    }

    std::string SAMFileWriter::getInLine ()
    {
        std::string str;
        while (!m_out_stream.eof())
        {
            std::getline(m_out_stream, str);
            return str;
        }

        return "";
    }

    long SAMFileWriter::getInStreamPosition () { return m_out_stream.tellg(); }
    long SAMFileWriter::getOutStreamPosition () { return m_out_stream.tellp(); }
    */

    void SAMFileWriter::printIosState ()
    {
        std::cout << "Is m_out_stream open? " << m_out_stream.is_open() << std::endl;
        std::cout << "good()=" << m_out_stream.good() << std::endl;
        std::cout << "eof()=" << m_out_stream.eof() << std::endl;
        std::cout << "fail()=" << m_out_stream.fail() << std::endl;
        std::cout << "bad()=" << m_out_stream.bad() << std::endl;
    }
}
