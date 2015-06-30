#include "core/file/ASCIIFileReader.h"
#include "core/file/ASCIIGZFileReader.h"
#include "VariantList.h"
#include "VCFFileReader.h"
#include "VCFVariant.h"

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

	VCFFileReader::VCFFileReader(const std::string& path)
	{
		static uint32_t s_vcf_id = 0; // An id that is set and auto increments when a new reader is created
		m_id = s_vcf_id;
		++s_vcf_id;
		setFileReader(path);
		Open();
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
	 * Later we might create a VCFHeader class
	 * but currently there is no need so for now
	 * this function just advances the pointer
	 * past the header.
	 */
	void VCFFileReader::readHeader()
	{
		std::string line;
		std::string headerEnd = "#CHROM";
		while (m_file_ptr->getNextLine(line) && strncmp(line.c_str(), headerEnd.c_str(), headerEnd.size()) != 0)
		{
		}
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
					variantPtrs.emplace_back(VCFVariant::BuildVariant(line, this->m_this_wk_ptr.lock(), this->m_vcf_parser));
				}
				if (regionPtr->getEndPosition() < linePosition) { break; } // if we have passed the end position of the region then stop looking for variants
			}
		}
		return variantPtrs;
	}

} // end namespace gwiz
