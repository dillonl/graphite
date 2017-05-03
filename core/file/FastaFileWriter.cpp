#include "FastaFileWriter.h"

#include <algorithm>

namespace graphite
{
    /*
	FastaWriter::FastaWriter(const std::string& header, const std::string& sequence) :
		m_header(header),
		m_sequence(sequence)
	{
	}

	FastaWriter::~FastaWriter()
	{
	}

	void FastaWriter::write(std::ostream& out)
	{
		if (this->m_header.size() > 0 && this->m_header[0] != '>')
		{
			out << ">";
		}
		out << this->m_header << std::endl;
		uint32_t lineLength = 80;
		uint32_t currentCharIndex = 0;
		for (uint32_t i = 0; i < this->m_sequence.size(); ++i)
		{
			if (i > 0 && (i % 80) == 0) { out << std::endl; }
			out.write(this->m_sequence.c_str() + i, 1);
		}
	}
    */

    FastaFileWriter::FastaFileWriter (std::vector< std::string >& headers, std::vector< std::string >& sequences) :
        m_headers(headers),
        m_sequences(sequences),
        m_opened(false)
    {
    }

    FastaFileWriter::~FastaFileWriter()
    {
    }

    bool FastaFileWriter::open (std::string fileName)
    {
        if (m_opened) { return true; }
        m_out_stream.open(fileName);
        m_opened = true;
    }

    void FastaFileWriter::write ()
    {
        uint32_t lineLength = 80;
        for (uint32_t i = 0; i < m_headers.size(); ++i)
        {
            for (uint32_t j = 0; j < m_headers[i].size(); ++j)
                m_out_stream.write(m_headers[i].c_str() + j, 1);

            m_out_stream << std::endl;

            for (uint32_t j = 0; j < m_sequences[i].size(); ++j)
            {
                if (j > 0 && (j % lineLength) == 0)
                    m_out_stream << std::endl;
                m_out_stream.write(m_sequences[i].c_str() + j, 1);
            }

            m_out_stream << std::endl;
        }
    }

    bool FastaFileWriter::close()
    {
        if (!m_opened) { return true; }
        m_out_stream.close();
        m_opened = false;

        return m_opened;
    }
}
