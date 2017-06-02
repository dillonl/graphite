#include "FastaFileWriter.h"

#include <algorithm>

namespace graphite
{
    //FastaFileWriter::FastaFileWriter (const std::string& path) :
    FastaFileWriter::FastaFileWriter () :
        //m_file_path(path),
        m_opened(false)
    {
    }

    FastaFileWriter::~FastaFileWriter()
    {
    }

    bool FastaFileWriter::open (std::string fileName)
    {
        if (m_opened) { return true; }
        //m_out_stream.open(m_file_path);
        m_out_stream.open(fileName, std::ios::app);
        m_opened = true;
    }

    void FastaFileWriter::write (const std::vector< std::string >& headers, const std::vector< std::string >& sequences)
    {
        uint32_t lineLength = 80;
        for (uint32_t i = 0; i < headers.size(); ++i)
        {
            for (uint32_t j = 0; j < headers[i].size(); ++j)
                m_out_stream.write(headers[i].c_str() + j, 1);

            m_out_stream << std::endl;

            for (uint32_t j = 0; j < sequences[i].size(); ++j)
            {
                if (j > 0 && (j % lineLength) == 0)
                    m_out_stream << std::endl;
                m_out_stream.write(sequences[i].c_str() + j, 1);
            }

            m_out_stream << std::endl;
        }
    }

    void FastaFileWriter::close()
    {
        if (!m_opened) { return; }
        m_out_stream.close();
        m_opened = false;
    }
}
