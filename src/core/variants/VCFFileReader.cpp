#include "VCFFileReader.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace gwiz
{
	VCFFileReader::VCFFileReader(const std::string& path, Region::SharedPtr region) :
		ASCIIFileReader(path)
	{
		Open(region);
	}

	VCFFileReader::VCFFileReader(const std::string& path) :
		ASCIIFileReader(path)
	{
		Open();
	}


	VCFFileReader::~VCFFileReader()
	{
	}

	void VCFFileReader::Open()
	{
		ASCIIFileReader::Open();
		this->m_current_position = this->m_file->const_data() + getHeaderSize(); // advance past the header
		this->m_end_position = this->m_file->const_data() + this->m_file_size; // set the end position to the eof
	}

	void VCFFileReader::Open(Region::SharedPtr region)
	{
		ASCIIFileReader::Open();
		if (region->getStartPosition() == 0 && region->getEndPosition() == 0)
		{
			region->setStartPosition(1);
			region->setEndPosition(std::numeric_limits< position >::max());
		}
		registerRegion(region);
		this->m_current_position = this->m_start_position;
	}

	void VCFFileReader::rewind()
	{
		this->m_current_position = this->m_start_position;
	}

	size_t VCFFileReader::getHeaderSize()
	{
		const char* start_pos = this->m_file->const_data();
		std::string end_header_value = "#CHROM";
		const char* end_pos = strstr(start_pos, end_header_value.c_str());
		end_pos = static_cast<const char*>(memchr(end_pos, '\n', std::numeric_limits< position >::max()));
		return (end_pos + 1) - start_pos;
	}

	/*
	 * Later we might create a VCFHeader class
	 * but currently there is no need so for now
	 * this function just advances the pointer
	 * past the header.
	 */
	void VCFFileReader::readHeader()
	{
		const char* line;
		std::string end_header_value = "#CHROM";
		do
		{
			line = getNextLine();
		} while (line != NULL && strncmp(line, end_header_value.c_str(), end_header_value.size()));
	}

	position VCFFileReader::getPositionFromLine(const char* line)
	{
		const char* posStart = static_cast<const char*>(memchr(line, '\t', 100)) + 1;
		const char* posEnd = static_cast<const char*>(memchr(posStart, '\t', 100));

		position pos;
		bool r = boost::spirit::qi::phrase_parse(
			posStart,
			posEnd,
			(
				boost::spirit::qi::uint_[boost::phoenix::ref(pos) = boost::spirit::qi::_1]
			),
			boost::spirit::ascii::space
			);
		return pos;
	}

	bool VCFFileReader::setRegionPositions(Region::SharedPtr regionPtr, const char* startLine, const char* endLine)
	{
		position startPositionTemp = std::numeric_limits< position >::max();
		position endPositionTemp = 0;
		const char* start = NULL;
		const char* end = NULL;
		std::string regionChrom = regionPtr->getReferenceID() + "\t";
		const char* regionChr = regionChrom.c_str();
		uint32_t regionChrSize = regionChrom.size();
		position regionStartPosition = regionPtr->getStartPosition();
		position regionEndPosition = regionPtr->getEndPosition();
		const char* startLineTemp = startLine;
		size_t fileSizeLeft = (this->m_end_file_ptr - startLineTemp);
		while (startLineTemp < endLine)
		{
			const char* endLineTemp = static_cast<const char*>(memchr(startLineTemp, '\n', fileSizeLeft));
			if (memcmp(startLineTemp, regionChr, regionChrSize) == 0)
			{
				position pos = getPositionFromLine(startLineTemp);
				bool withinStartPosition = pos >= regionStartPosition;
				if (withinStartPosition && pos < startPositionTemp) // check start position
				{
					startPositionTemp = pos;
					start = startLineTemp;
				}
				if (withinStartPosition && pos <= regionEndPosition && pos > endPositionTemp) // check end position
				{
					endPositionTemp = pos;
					end = startLineTemp;
				}
			}
			startLineTemp = endLineTemp + 1;
		}

		this->m_region_mutex.lock();
		if (startPositionTemp != std::numeric_limits< position >::max())
		{
			this->m_start_position = (this->m_start_position != NULL && getPositionFromLine(this->m_start_position) < startPositionTemp) ? this->m_start_position : start;
		}
		if (endPositionTemp != 0)
		{
			this->m_end_position = (this->m_end_position != NULL && getPositionFromLine(this->m_end_position) > endPositionTemp) ? this->m_end_position : end;
			this->m_end_position = static_cast< const char* >(memchr(this->m_end_position, '\n', std::numeric_limits< size_t >::max())) + 1;
		}
		this->m_region_mutex.unlock();
		return true;
	}

	bool VCFFileReader::registerRegion(Region::SharedPtr region)
	{
		size_t headerSize = getHeaderSize();
		const char* startPosition = this->m_file->const_data() + headerSize;
		const char* endPosition = this->m_end_file_ptr;

		size_t distance = endPosition - startPosition;
		uint32_t ONE_HUNDRED_MEGA_BYTES = 100000000;
		uint32_t maxNumberOfThreads = std::max< uint32_t >(std::thread::hardware_concurrency(), 10);
		uint32_t numberOfThreads = std::min<uint32_t>(maxNumberOfThreads, std::max<uint32_t>((this->m_file_size/ONE_HUNDRED_MEGA_BYTES), 1)); // 1 thread per 100mb of data (min = 1 max = maxNumberofthreads)
		size_t chunk = (distance/numberOfThreads); // size of the chunk of memory to work on
		std::vector< std::thread > threads;
		const char* startLine = startPosition;
		const char* endLine;
		for (uint32_t i = 0; i < numberOfThreads; ++i)
		{
			size_t tempFileSize = endPosition - startLine + chunk;
			endLine = (i == numberOfThreads - 1) ? endPosition : static_cast<const char*>(memchr(startLine + chunk, '\n', tempFileSize));
			threads.push_back(std::thread(&VCFFileReader::setRegionPositions, this, region, startLine, endLine));
			startLine = endLine + 1;
		}
		for (auto iter = threads.begin(); iter != threads.end(); ++iter)
		{
			(*iter).join();
		}
		if (this->m_start_position == NULL || this->m_end_position == NULL)
		{
			throw "Region " + region->getRegionString() + " was not found.";
		}
		return true;
	}

	void VCFFileReader::printToVCF(std::ostream& out)
	{
		throw "VCFFileReader printToVCF not implemented";
	}

} // end namespace gwiz
