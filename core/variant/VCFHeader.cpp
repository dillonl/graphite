#include "VCFHeader.h"

namespace graphite
{
	VCFHeader::VCFHeader()
	{
		m_header_lines.emplace_back("##fileformat=VCFv4.2");
		m_header_lines.emplace_back("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
		m_header_lines.emplace_back("##INFO=<ID=TC,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">");
		m_header_lines.emplace_back("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth including non-counted low quality reads\">");
		m_header_lines.emplace_back("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		m_header_lines.emplace_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE");
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
		return header;
	}
}
