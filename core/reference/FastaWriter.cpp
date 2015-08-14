#include "FastaWriter.h"

#include <fstream>
#include <algorithm>

namespace graphite
{
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
}
