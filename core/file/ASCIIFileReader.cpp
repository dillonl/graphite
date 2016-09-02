#include "ASCIIFileReader.h"

#include <memory>
#include <stdexcept>
#include <algorithm>  // for std::find
#include <cstring>

namespace graphite
{
	ASCIIFileReader::ASCIIFileReader(const std::string& path) :
		IFile(path),
		m_start_position(NULL),
		m_end_position(NULL)
		//m_end_file_ptr(0)
	{
		// this->m_line.reserve(10000);
	}

	ASCIIFileReader::~ASCIIFileReader()
	{
		if (this->m_opened) // if open close when deleted
		{
			this->Close();
		}
	}

	void ASCIIFileReader::Open()
	{
		if (this->m_opened) { return; } // only open one time
		m_file = std::make_shared< std::ifstream >(m_file_path, std::fstream::in);
		// if (!this->fileExists(this->m_file_path))
		if (!m_file->good())
		{
            throw std::ios_base::failure("File not found: " + this->m_file_path);
		}
		// this->m_file = std::make_shared<boost::iostreams::mapped_file>(this->m_file_path.c_str(), boost::iostreams::mapped_file::readonly);
		// this->m_file_size = this->m_file->size();
		// this->m_current_position = this->m_file->const_data();
		// this->m_end_file_ptr = (this->m_current_position + this->m_file_size);
		this->m_opened = true;
		// this->m_end_position = this->m_file->const_data() + this->m_file_size; // set the end position to the eof
		// this->m_end_position = this->m_end_file_ptr;
	}

	void ASCIIFileReader::Close()
	{
		if (this->m_opened)
		{
			this->m_file->close();
			this->m_opened = false;
		}
	}
}
