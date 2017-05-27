#include "BamFileReader.h"

namespace graphite
{
    BamFileReader::BamFileReader (const std::string& bamPath) :
        m_bam_path(bamPath),
        m_opened(false)
    {
    }

    BamFileReader::~BamFileReader ()
    {
    }

    void BamFileReader::open ()
    {
        if (m_opened) { return; }
        m_bam_reader.Open(m_bam_path);
        m_opened = true;
    }

    std::string BamFileReader::getSamHeader ()
    {
        return m_bam_reader.GetHeaderText();
    }

    void BamFileReader::close ()
    {
        if (m_opened)
        {
            m_bam_reader.Close();
            m_opened = false;
        }
    }
}
