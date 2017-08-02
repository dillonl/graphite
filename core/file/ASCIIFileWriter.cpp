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
			close();
		}
	}

	bool ASCIIFileWriter::open()
	{
		if (m_opened) { return false; }
		this->m_out_stream.open(this->m_file_path);
		m_opened = true;
		return true;
	}

	void ASCIIFileWriter::close()
	{
		if (!m_opened) { return; }
		this->m_out_stream.close();
		m_opened = false;
	}

	bool ASCIIFileWriter::write(const char* data, size_t dataLength)
	{
		if (!m_opened) { return false; }
		this->m_out_stream.write(data, dataLength);
        m_out_stream << std::endl;
		m_opened = true;
	}
}
