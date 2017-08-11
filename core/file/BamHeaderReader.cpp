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

    void BamHeaderReader::addPathHeaderToSamHeader (std::string graphPathHeader, int graphPathLength)
    {
        BamTools::SamSequence samSequence = BamTools::SamSequence(graphPathHeader, graphPathLength);
        m_sam_header.Sequences.Add(samSequence);
    }

    void BamHeaderReader::addReadGroupsToSamHeader ()
    {
        m_sam_header.ReadGroups.Add("REF");
        m_sam_header.ReadGroups.Add("ALT");
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

    bool BamHeaderReader::containsReadGroup (std::string readGroup)
    {
        return m_sam_header.ReadGroups.Contains(readGroup);
    }
}
