#include "core/file/ASCIIFileReader.h"
#include "core/file/ASCIIGZFileReader.h"
#include "core/region/Region.h"
#include "VariantList.h"
#include "VCFFileReader.h"
#include "Variant.h"

#include <fstream>
#include <iostream>
#include <sstream>

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
		while (m_file_ptr->getNextLine(line) && line.size() > 0 && line.c_str()[0] == '#')
		{
			this->m_vcf_header->addHeaderLine(line);
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
		/*
		const char* posStart = static_cast<const char*>(memchr(line, '\t', 100)) + 1;
		const char* posEnd = static_cast<const char*>(memchr(posStart, '\t', 100));

		while (posStart != posEnd)
		{

		}
		*/

		/*
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
		*/
	}

	std::vector< IVariant::SharedPtr > VCFFileReader::getVariantsInRegion(Region::SharedPtr regionPtr)
	{
		std::vector< IVariant::SharedPtr > variantPtrs;
		std::string regionReferenceIDWithTab = regionPtr->getReferenceID() + "\t";
		std::string line;
		uint32_t count = 0;
		while (this->m_file_ptr->getNextLine(line))
		{
			if (memcmp(regionReferenceIDWithTab.c_str(), line.c_str(), regionReferenceIDWithTab.size()) == 0) // if we are in the correct reference (chrom)
			{
				position linePosition = getPositionFromLine(line.c_str());
				if ((regionPtr->getStartPosition() <= linePosition && linePosition <= regionPtr->getEndPosition()))
				{
					variantPtrs.emplace_back(Variant::BuildVariant(line, this->m_reference_ptr, m_max_allowed_allele_size));
				}
				if (regionPtr->getEndPosition() < linePosition) { break; } // if we have passed the end position of the region then stop looking for variants
			}
		}
		return variantPtrs;
	}

} // end namespace graphite
