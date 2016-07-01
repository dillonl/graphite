#include "VCFHeader.h"
#include "core/alignment/Sample.hpp"

#include <unordered_set>

namespace graphite
{
	VCFHeader::VCFHeader()
	{
		m_header_lines.emplace_back("##fileformat=VCFv4.2");
		m_header_lines.emplace_back("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
		m_header_lines.emplace_back("##INFO=<ID=DP4,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">");
	}

	VCFHeader::~VCFHeader() {}

	void VCFHeader::addHeaderLine(const std::string& headerLine)
	{
		m_header_lines.emplace_back(headerLine);
	}

	std::string VCFHeader::getHeader()
	{
		std::string header = "";
		for (std::string& headerLine : m_header_lines)
		{
			header += headerLine + "\n";
		}
		std::string sampleHeaderString;
		std::unordered_set< std::string > sampleNames;
		for (auto samplePtr : m_sample_ptrs)
		{
			auto iter = sampleNames.find(samplePtr->getName());
			if (iter == sampleNames.end())
			{
				sampleHeaderString += "\t" + samplePtr->getName();
				sampleNames.emplace(samplePtr->getName());
			}
		}
		header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + sampleHeaderString + "\n";
		return header;
	}

	void VCFHeader::registerSample(std::shared_ptr< Sample > samplePtr)
	{
		m_sample_ptrs.emplace_back(samplePtr);
	}
}
