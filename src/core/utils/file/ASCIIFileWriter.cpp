#include "ASCIIFileWriter.h"

namespace gwiz
{
	ASCIIFileWriter::ASCIIFileWriter(const std::string& path) :
		IFile(path)
	{
	}

	ASCIIFileWriter::~ASCIIFileWriter()
	{
	}

	void ASCIIFileWriter::Write(const char* output, const size_t size)
	{
		this->m_out_stream.write(output, size);
	}

	void ASCIIFileWriter::Open()
	{
		this->m_out_stream.open(this->m_file_path);
	}

	void ASCIIFileWriter::Close()
	{
		this->m_out_stream.close();
	}
}
