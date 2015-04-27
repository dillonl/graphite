#include "ASCIIGZFileReader.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

namespace gwiz
{
	ASCIIGZFileReader::ASCIIGZFileReader(const std::string& path) :
		IFile(path)
	{
	}

	ASCIIGZFileReader::~ASCIIGZFileReader()
	{
	}

	void ASCIIGZFileReader::Open()
	{
		m_ifstream_ptr = std::make_shared< std::ifstream >(this->m_file_path, std::ios_base::in | std::ios_base::binary);
		m_in_stream_ptr = std::make_shared< boost::iostreams::filtering_istream >();
		m_in_stream_ptr->push(boost::iostreams::gzip_decompressor());
		m_in_stream_ptr->push(*m_ifstream_ptr);
	}

	void ASCIIGZFileReader::Close()
	{
		m_ifstream_ptr->close();
	}
}
