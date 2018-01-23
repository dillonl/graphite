#include "ASCIIFileWriter.h"

namespace graphite
{
	ASCIIFileWriter::ASCIIFileWriter(const std::string& path) :
		IFileWriter(FileType::ASCII),
		m_file_path(path),
		m_opened(false)
	{
	}

	ASCIIFileWriter::~ASCIIFileWriter()
	{
		if (m_opened)
		{
			std::cout << "close ascii file writer" << std::endl;
			close();
		}
	}

	bool ASCIIFileWriter::write(const char* data, size_t dataLength)
	{
		if (m_opened)
		{
			this->m_out_stream.write(data, dataLength);
			m_opened = true;
			return true;
		}
		else
		{
			return false;
		}
	}

	bool ASCIIFileWriter::open()
	{
		if (!m_opened)
		{
			this->m_out_stream.open(this->m_file_path);
			m_opened = true;
			return true;
		}
		else
		{
			return false;
		}
	}

	void ASCIIFileWriter::close()
	{
		if (m_opened)
		{
			this->m_out_stream.close();
			m_opened = false;
		}
	}
}
