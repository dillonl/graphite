#include "BamHeaderReader.h"

namespace graphite
{
    BamHeaderReader::BamHeaderReader (const std::string& bamPath) :
        m_bam_path(bamPath),
        m_opened(false)
    {
    }

    BamHeaderReader::~BamHeaderReader ()
    {
    }

    void BamHeaderReader::open ()
    {
        if (m_opened) { return; }
        m_bam_reader.Open(m_bam_path);
        m_opened = true;
        m_sam_header = m_bam_reader.GetHeader();
    }

    void BamHeaderReader::createPathHeaderVector (std::vector< std::string > graphPathHeaders, std::vector< int > graphPathLengths)
    {
        // Make sure the header and length vectors are the same size.
        assert(graphPathHeaders.size() == graphPathLengths.size());
        for (int i = 0; i < graphPathHeaders.size(); ++i)
        {
            BamTools::SamSequence samSequence = BamTools::SamSequence(graphPathHeaders[i], graphPathLengths[i]);
            m_sam_sequences.push_back(samSequence);
        }
    }

    void BamHeaderReader::addPathHeadersToSamHeader ()
    {
        m_sam_header.Sequences.Add(m_sam_sequences);
    }

    void BamHeaderReader::close ()
    {
        if (m_opened)
        {
            m_bam_reader.Close();
            m_opened = false;
        }
    }

    std::string BamHeaderReader::getOriginalSamHeader ()
    {
        return m_bam_reader.GetHeaderText();
    }
    
    std::string BamHeaderReader::getModifiedSamHeader ()
    {
        return m_sam_header.ToString();
    }
}
