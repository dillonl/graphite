#ifndef GRAPHITE_ASCIIGZFILEREADER_H
#define GRAPHITE_ASCIIGZFILEREADER_H

#include "IFile.h"

#include <memory>
#include <iostream>
#include <fstream>

/*
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/seek.hpp>
*/
/* #include <boost/iostreams/positioning.hpp> */
namespace graphite
{
	class ASCIIGZFileReader : public IFile
	{
	public:
		typedef std::shared_ptr< ASCIIGZFileReader > SharedPtr;
		ASCIIGZFileReader(const std::string& path);
		~ASCIIGZFileReader();

		virtual void Open() override;

		void Close() override;

		inline bool getNextLine(std::string& line) override
		{
			/*
			if (!this->m_opened) { return false; }
			auto value = (bool)std::getline(*this->m_in_stream_ptr, line);
			return value;
			*/
			return false;
		}

	private:
		std::shared_ptr< std::iostream > m_iostream_ptr;
		std::shared_ptr< std::ifstream > m_ifstream_ptr;
		/* std::shared_ptr< boost::iostreams::filtering_istream > m_in_stream_ptr; */
	};
}

#endif //GRAPHITE_ASCIIGZFILEREADER_H
