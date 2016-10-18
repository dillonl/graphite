#include "core/file/ASCIIFileReader.h"
#include "core/file/ASCIIGZFileReader.h"
#include "core/util/ThreadPool.hpp"
#include "core/region/Region.h"
#include "VariantList.h"
#include "VCFFileReader.h"
#include "Variant.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <deque>
#include <future>
#include <atomic>
#include <memory>

namespace graphite
{

	VCFFileReader::VCFFileReader(const std::string& path, IReference::SharedPtr referencePtr, uint32_t maxAlleleSize) :
		m_path(path),
		m_reference_ptr(referencePtr),
		m_max_allowed_allele_size(maxAlleleSize)
	{
		m_vcf_header = std::make_shared< VCFHeader >();
		static uint32_t s_vcf_id = 0; // An id that is set and auto increments when a new reader is created
		m_id = s_vcf_id;
		++s_vcf_id;
		setFileReader(m_path);
		Open();
	}

	VCFFileReader::VCFFileReader(const std::string& path) :
		m_path(path)
	{
		m_vcf_header = std::make_shared< VCFHeader >();
		setFileReader(m_path);
		Open();
	}

	std::string VCFFileReader::getFilePath()
	{
		return this->m_path;
	}

	void VCFFileReader::setFileReader(const std::string& path)
	{
		std::string fileExtension = path.substr(path.find_last_of('.') + 1);
		if (strcmp(fileExtension.c_str(), "gz") == 0) // if it's a gz file
		{
			m_file_ptr = std::make_shared< ASCIIGZFileReader >(path);
		}
		else // if it's not gz'd
		{
			m_file_ptr = std::make_shared< ASCIIFileReader >(path);
		}
	}


	VCFFileReader::~VCFFileReader()
	{
	}

	void VCFFileReader::Open()
	{
		m_file_ptr->Open();
		readHeader();
	}

	/*
	 * Set the reader's vcfheader with this file's vcf header.
	 */
	void VCFFileReader::readHeader()
	{
		std::string line;
		std::string headerEnd = "#CHROM";
		while (m_file_ptr->getNextLine(line) && line.size() > 0)
		{
			this->m_vcf_header->addHeaderLine(line);
			if (strncmp(headerEnd.c_str(), line.c_str(), headerEnd.size()) == 0) { break; }
		}
	}

	std::vector< Region::SharedPtr > VCFFileReader::GetAllRegionsInVCF(const std::string& vcfPath)
	{
		std::vector< Region::SharedPtr > regionPtrs;
		std::string line;
		std::string currentRegion = "";
		auto vcfFileReaderPtr = std::make_shared< VCFFileReader >(vcfPath);
		while (vcfFileReaderPtr->m_file_ptr->getNextLine(line))
		{
			auto region = line.substr(0, line.find("\t"));
			if (currentRegion.compare(region) != 0)
			{
				auto regionPtr = std::make_shared< Region >(region);
				regionPtrs.emplace_back(regionPtr);
				currentRegion = region;
			}
		}
		return regionPtrs;
	}

	position VCFFileReader::getPositionFromLine(const char* line)
	{
		const char* tmpLine = line;
		size_t posSize = 0;
		while (*tmpLine != '\n')
		{
			if (*tmpLine == '\t')
			{
				break;
			}
			++tmpLine;
		}
		char* endPtr;
		auto pos = strtol(tmpLine, &endPtr, 10);
		return pos;
	}

	std::vector< IVariant::SharedPtr > VCFFileReader::getVariantsInRegion(Region::SharedPtr regionPtr)
	{
		std::vector< IVariant::SharedPtr > variantPtrs;
		std::string regionReferenceIDWithTab = regionPtr->getReferenceID() + "\t";
		std::string line;
		uint32_t count = 0;
		bool wasInRegion = false;
		// this->m_file_ptr->setFilePosition(findRegionStartPosition(regionPtr));
		// std::cout << "region not yet found" << std::endl;
		while (this->m_file_ptr->getNextLine(line))
		{
			if (memcmp(regionReferenceIDWithTab.c_str(), line.c_str(), regionReferenceIDWithTab.size()) == 0) // if we are in the correct reference (chrom)
			{
				position linePosition = getPositionFromLine(line.c_str());
				if ((regionPtr->getStartPosition() <= linePosition && linePosition <= regionPtr->getEndPosition()))
				{
					variantPtrs.emplace_back(Variant::BuildVariant(line, this->m_reference_ptr, m_max_allowed_allele_size));
					wasInRegion = true;
					continue;
				}
				if (regionPtr->getEndPosition() < linePosition) { break; } // if we have passed the end position of the region then stop looking for variants
			}
			if (wasInRegion)
			{
				break;
			}
		}
		return variantPtrs;
	}

	uint64_t VCFFileReader::findRegionStartPosition(Region::SharedPtr regionPtr)
	{
		std::ifstream in(this->m_path, std::ifstream::ate | std::ifstream::binary);
		auto fileSize = in.tellg();
		in.close();
		std::deque< std::shared_ptr< std::future< uint64_t > > > futureFunctions;
		std::shared_ptr< std::atomic< bool > > posFound = std::make_shared< std::atomic< bool > >(false);
		uint32_t numThreads = 20;
		uint32_t partSize = fileSize / numThreads;
		uint64_t seekPosition = 0;
		for (auto i = 0; i < numThreads; ++i)
		{
			// std::cout << "starting thread: " << i << " at seek position: " << seekPosition << std::endl;
			auto funct = std::bind(&VCFFileReader::getPositionFromFile, this, seekPosition, seekPosition + partSize, posFound, regionPtr, this->m_path);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			seekPosition += partSize;
		}
		uint64_t pos = 0;
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
			{
				auto tmpPos = futureFunct->get();
				if (tmpPos > 0)
				{
					pos = tmpPos;
					break;
				}
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}

		return pos;
	}

	uint64_t VCFFileReader::getPositionFromFile(uint64_t seekPosition, uint64_t endSeekPosition, std::shared_ptr< std::atomic< bool > > posFound, Region::SharedPtr regionPtr, std::string path)
	{
		static std::mutex l;

		std::ifstream f;
		f.open(path.c_str(), std::ifstream::in);

		f.seekg(endSeekPosition, std::ifstream::beg);
		std::string line;

		bool newChrom = false;

		if (seekPosition > 0)
		{
			// rewind by one line
			while (f.unget() && f.peek() != '\n')
			{
				std::lock_guard< std::mutex > lock(l);
				char tmp = f.peek();
				// std::cout << tmp << std::endl;
			}
		}

		std::getline(f, line); // fastforward to the eol
		endSeekPosition = f.tellg();
		f.seekg(seekPosition, std::ifstream::beg);

		std::string regionReferenceIDWithTab = regionPtr->getReferenceID() + "\t";
		position startPosition = regionPtr->getStartPosition();
		uint64_t returnFilePos = 0;
		uint64_t tmpFilePos = f.tellg();
		std::getline(f, line); // get the partial line
		uint64_t counter = 0;
		while (std::getline(f, line) && f.tellg() < endSeekPosition)
		{
			{
				std::lock_guard< std::mutex > lock(l);
				// std::cout << line << std::endl;
			}
			if (line.size() > 0 && line[0] == '#') { continue; }
			if (counter++ % 10 == 0 && posFound->load()) { break; } // another thread located the position before this thread
			if (memcmp(regionReferenceIDWithTab.c_str(), line.c_str(), regionReferenceIDWithTab.size()) == 0) // if we are in the correct reference (chrom)
			{
				auto currentPosition = VCFFileReader::getPositionFromLine(line.c_str());
				if (currentPosition >= startPosition)
				{
					returnFilePos = tmpFilePos;
					f.seekg(returnFilePos, std::ifstream::beg);
					std::getline(f, line); // get the partial line
					{
						std::lock_guard< std::mutex > lock(l);
						// std::cout << "found position: " << currentPosition << " " << returnFilePos << std::endl;
						// std::cout << "found line: " << line << std::endl;
					}
					*posFound = true;
					break;
				}
			}
			tmpFilePos = f.tellg();
		}
		f.close();
		{
			std::lock_guard< std::mutex > lock(l);
			// std::cout << "finished: " << returnFilePos << std::endl;
		}
		return returnFilePos;
	}

} // end namespace graphite
