#include "ASCIIGZFileReader.h"

#include <zconf.h>


namespace graphite
{
	ASCIIGZFileReader::ASCIIGZFileReader(const std::string& path) :
		IFile(path)
	{
		Open();
	}

	ASCIIGZFileReader::~ASCIIGZFileReader()
	{
		Close();
	}

	void ASCIIGZFileReader::Open()
	{
		if (m_opened) { return; }
		if (!this->fileExists(this->m_file_path))
		{
            throw std::ios_base::failure("File not found: " + this->m_file_path);
		}
		m_opened = true;
		m_in_stream_ptr = std::make_shared< igzstream >();
		m_in_stream_ptr->open(this->m_file_path.c_str());
		/*
		m_ifstream_ptr = std::make_shared< std::ifstream >(this->m_file_path, std::ios_base::in | std::ios_base::binary);
		m_in_stream_ptr = std::make_shared< boost::iostreams::filtering_istream >();
		m_in_stream_ptr->push(boost::iostreams::gzip_decompressor());
		m_in_stream_ptr->push(*m_ifstream_ptr);
		*/
	}

	void ASCIIGZFileReader::Close()
	{
		if (!m_opened) { return; }
		m_in_stream_ptr->close();
		m_opened = false;
	}
}
