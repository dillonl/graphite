#include "GZASCIIFileReader.h"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <memory>
#include <stdexcept>
#include <algorithm>  // for std::find
#include <cstring>

namespace gwiz
{
	GZASCIIFileReader::GZASCIIFileReader(std::string& path) :
		IFile(path),
		m_start_position(NULL),
		m_end_position(NULL)
		//m_end_file_ptr(0)
	{
		// this->m_line.reserve(10000);
	}

	GZASCIIFileReader::~GZASCIIFileReader()
	{
		if (this->m_opened) // if open close when deleted
		{
			this->Close();
		}
	}

	void GZASCIIFileReader::Open()
	{
		if (this->m_opened) { return; } // only open one time
		if (!this->fileExists(this->m_file_path))
		{
            throw std::ios_base::failure("File not found: " + this->m_file_path);
		}

this->m_file = std::make_shared< boost::iostreams::mapped_file_source >(this->m_file_path.c_str());

		this->m_file = std::make_shared< boost::iostreams::mapped_file >(this->m_file_path.c_str(), boost::iostreams::mapped_file::readonly);
        this->m_text_stream = std::make_shared< boost::iostreams::stream< io::mapped_file_source > >(this->m_file);
        this->m_filtering_stream = std::make_shared< boost::iostreams::filtering_istream > fs;
        this->m_filtering_stream->push(boost::iostreams::gzip_decompressor{});
        this->m_filtering_stream->push(textstream);

		this->m_file_size = this->m_file->size();
		this->m_current_position = this->m_file->const_data();
		this->m_end_file_ptr = (this->m_current_position + this->m_file_size);
		this->m_opened = true;
		// this->m_end_position = this->m_end_file_ptr;
	}

	void GZASCIIFileReader::Close()
	{
		if (this->m_opened)
		{
			this->m_file->close();
			this->m_opened = false;
		}
	}
}
