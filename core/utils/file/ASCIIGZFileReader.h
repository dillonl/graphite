#ifndef GWIZ_ASCIIGZFILEREADER_H
#define GWIZ_ASCIIGZFILEREADER_H

#include "IFile.h"

#include <memory>
#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/seek.hpp>
/* #include <boost/iostreams/positioning.hpp> */
namespace gwiz
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
			auto value = (bool)std::getline(*this->m_in_stream_ptr, line);
			line += "\n";
			return value;
		}

		/*
		inline const char* getNextLine() override
		{
			std::cout << "hello" << std::endl;
			for(std::string str; std::getline(*this->m_iostream_ptr, str); )
			{
				std::cout << "Processed line " << str << '\n';
			}
			return nullptr;
		}
		*/

	private:
		std::shared_ptr< std::iostream > m_iostream_ptr;
		std::shared_ptr< std::ifstream > m_ifstream_ptr;
		std::shared_ptr< boost::iostreams::filtering_istream > m_in_stream_ptr;
	};
}

#endif //GWIZ_ASCIIGZFILEREADER_H
