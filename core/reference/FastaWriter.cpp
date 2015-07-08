#include "FastaWriter.h"

#include <fstream>
#include <algorithm>

namespace gwiz
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
		while (currentCharIndex < this->m_sequence.size())
		{
			out.write((this->m_sequence.c_str() + currentCharIndex), lineLength);
			currentCharIndex += lineLength;
			lineLength = std::min< uint32_t >(lineLength, std::max< uint32_t >((uint32_t)(this->m_sequence.size() - currentCharIndex), 0));
		}
	}
}
