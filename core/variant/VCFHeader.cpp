#include "VCFHeader.h"
#include "core/alignment/Sample.hpp"

#include <unordered_set>

namespace graphite
{
	VCFHeader::VCFHeader()
	{
		// m_header_lines.emplace_back("##fileformat=VCFv4.2");
		// m_header_lines.emplace_back("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
		// m_header_lines.emplace_back("##INFO=<ID=DP4,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">");
	}

	VCFHeader::~VCFHeader() {}

	void VCFHeader::addHeaderLine(const std::string& headerLine)
	{
		m_header_lines.emplace_back(headerLine);
	}

	std::string VCFHeader::getHeader()
	{
		std::string sampleHeaderString;
		std::unordered_set< std::string > sampleNames;
		for (auto samplePtr : m_sample_ptrs)
		{
			auto iter = sampleNames.find(samplePtr->getName());
			if (iter == sampleNames.end())
			{
				sampleHeaderString += "\t" + samplePtr->getName();
				sampleNames.emplace(samplePtr->getName()); // this is to make sure we don't add duplicate samples to the header
			}
		}

		std::string header = "";
		std::string firstLinePrefix = "##fileformat";
		std::string headerEnd = "#CHROM";
		std::string headerInfo = "##INFO";
		std::string headerFilter = "##FILTER";
		std::string headerFormat = "##FORMAT";
		bool insertedInfo = false;
		bool firstLine = true;
		for (std::string& headerLine : m_header_lines)
		{
			if (firstLine && headerLine.compare(0, firstLinePrefix.size(), headerLine) != 0) // the first line should alway be this. If it's not then add it
			{
				header += "##fileformat=VCFv4.2";
			}
			if (!insertedInfo &&
				(headerLine.compare(0, headerInfo.size(), headerInfo) == 0 ||
				 headerLine.compare(0, headerFilter.size(), headerFilter) == 0 ||
				 headerLine.compare(0, headerFormat.size(), headerFormat) == 0 ||
				 headerLine.compare(0, headerEnd.size(), headerEnd) == 0)) // we need to put in the graphite info fields
			{
				header += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n";
				header += "##INFO=<ID=DP4,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">\n";
				insertedInfo = true;
			}
			else if (headerLine.compare(0, headerEnd.size(), headerEnd) == 0) // if we are at the final point of the vcf header then insert the sample names
			{
				header += headerLine + sampleHeaderString + "\n";
			}
			else
			{
				header += headerLine + "\n";
			}
			firstLine = false;
		}
		return header;
	}

	void VCFHeader::registerReferencePath(const std::string& referencePath)
	{
		m_reference_path = referencePath;
	}

	void VCFHeader::registerSample(std::shared_ptr< Sample > samplePtr)
	{
		m_sample_ptrs.emplace_back(samplePtr);
	}
}
