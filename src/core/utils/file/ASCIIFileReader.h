#ifndef GWIZ_ASCIIFILEREADER_H
#define GWIZ_ASCIIFILEREADER_H

#include <fcntl.h>    /* For O_RDWR */
#include <unistd.h>   /* For open(), creat() */
#include <map>

#include <iostream>   // for std::cout
#include <boost/iostreams/device/mapped_file.hpp> // for mmap

#include "IFile.h"

namespace gwiz
{

	class ASCIIFileReader : public IFile
	{
	public:
		typedef std::shared_ptr< ASCIIFileReader > SharedPtr;
		ASCIIFileReader(const std::string& path);
		~ASCIIFileReader() override;

		/*
		 * Opens (and sets the size of the file) the file using
		 * the constructor's path. If the file does not exist
		 * then an exception is thrown.
		 */
		virtual void Open() override;
		/*
		 * Closes the opened file. If the file is not
		 * opened then this function does nothing.
		 */
		void Close() override;

		/*
		 * Returns a handle to the m_line pointer.
		 * This may need to be modified to support
		 * multithreaded functionality.
		 * Advances file position.
		 */
		inline const char* getNextLine()
		{
			if (!this->m_opened || this->m_current_position >= this->m_end_position) { return NULL; }
			const char* current_line = this->m_current_position;
			const char* line = static_cast< const char* >(memchr(this->m_current_position, '\n', (this->m_end_file_ptr - this->m_current_position)));
			size_t lineSize = (line - this->m_current_position) + 1;
			this->m_current_position += lineSize;
			return current_line;
		}

	protected:
		std::shared_ptr<boost::iostreams::mapped_file> m_file;
		const char* m_start_position;
		const char* m_end_position;
		const char* m_current_position;

		const char* m_end_file_ptr;
	};

} // end namespace gwiz

#endif //GWIZ_ASCIIFILEREADER_H
