#ifndef GWIZ_GZASCIIFILEREADER_H
#define GWIZ_GZASCIIFILEREADER_H

#include <fcntl.h>    /* For O_RDWR */
#include <unistd.h>   /* For open(), creat() */
#include <map>

#include <iostream>   // for std::cout
#include <boost/iostreams/device/mapped_file.hpp> // for mmap

#include "IFile.h"

namespace gwiz
{

	class GZASCIIFileReader : public IFile
	{
	public:
		typedef std::shared_ptr<GZASCIIFileReader> SharedPtr;
		GZASCIIFileReader(std::string& path);
		~GZASCIIFileReader() override;

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
		inline uintmax_t CountLines(std::string& fname)
		{
			static const auto BUFFER_SIZE = 16*1024;
			int fd = open(fname.c_str(), O_RDONLY);
			if(fd == -1)

			posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL

			char buf[BUFFER_SIZE + 1];
			uintmax_t lines = 0;

			while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
			{
				if(bytes_read == (size_t)-1)
					//handle_error("read failed");
				if (!bytes_read)
					break;

				for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
					++lines;
			}

			return lines;
		}

		inline uintmax_t CountLines2(std::string& fname)
		{
			boost::iostreams::mapped_file mmap(fname.c_str(), boost::iostreams::mapped_file::readonly);
			auto f = mmap.const_data();
			auto l = f + mmap.size();

			uintmax_t m_numLines = 0;
			while (f && f!=l)
				if ((f = static_cast<const char*>(memchr(f, '\n', l-f))))
					m_numLines++, f++;

			return m_numLines;
		}

		inline void setLinePositions()
		{
			std::cout << "starting" << std::endl;
			if (this->m_opened || (this->m_end_file_ptr <= this->m_data_ptr)) { return; }
			this->m_line_positions.clear();
			const char* line = this->m_data_ptr;
			do
			{
				line = static_cast<const char*>(rawmemchr(line, '\n'));
				m_line_positions.push_back(line);
				++line; // move one past the \n
			} while (line < this->m_end_file_ptr);
			std::cout << "ending: " << m_line_positions.size() * sizeof(char*) << std::endl;
		}
		*/

		/*
		 * Returns a handle to the m_line pointer.
		 * This may need to be modified to support
		 * multithreaded functionality.
		 */
		inline const char* getNextLine()
		{
			if (!this->m_opened || this->m_current_position > this->m_end_position) { return NULL; }
			const char* current_line = this->m_current_position;
			const char* line = static_cast< const char* >(memchr(this->m_current_position, '\n', (this->m_end_file_ptr - this->m_current_position)));
			size_t lineSize = (line - this->m_current_position) + 1;
			this->m_current_position += lineSize;
			return current_line;
		}

	protected:
		std::shared_ptr< boost::iostreams::mapped_file_source > m_file;
		std::shared_ptr< boost::iostreams::stream< io::mapped_file_source > > m_text_stream;
        std::shared_ptr< boost::iostreams::filtering_istream > m_filtering_stream;
		const char* m_start_position;
		const char* m_end_position;
		const char* m_current_position;

		const char* m_end_file_ptr;
	};

} // end namespace gwiz

#endif //GWIZ_ASCIIFILEREADER_H
