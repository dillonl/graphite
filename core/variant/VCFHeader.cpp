#include "VCFHeader.h"
#include "core/util/Utility.h"

#include <unordered_set>
#include <algorithm>
#include <string.h>

namespace graphite
{
	VCFHeader::VCFHeader(const std::vector< std::string >& vcfHeaderLines)
	{
		std::string headerEnd = "#CHROM";
		std::string formatString = "##FORMAT";
		bool addedFormat = false;
		for (auto headerLine : vcfHeaderLines)
		{
			std::string upperLine = headerLine;
            std::transform(upperLine.begin(), upperLine.end(),upperLine.begin(), ::toupper);
			if (strncmp(headerEnd.c_str(), upperLine.c_str(), headerEnd.size()) == 0)
			{
				if (!addedFormat)
				{
					addedFormat = true;
					addFormatToHeader();
				}
				setColumns(headerLine);
			}
			else if (!addedFormat && strncmp(formatString.c_str(), upperLine.c_str(), formatString.size()) == 0)
			{
				addedFormat = true;
				addFormatToHeader();
				m_lines.emplace_back(headerLine);
			}
			else
			{
				m_lines.emplace_back(headerLine);
			}
		}
	}

	VCFHeader::~VCFHeader() {}

	void VCFHeader::setColumns(const std::string& headerString)
	{
		std::vector< std::string > requiredColumns = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};

		m_columns.clear();
		std::vector< std::string > headerSplit;
		split(headerString, '\t', headerSplit);
		bool invalidHeader = false;
		uint32_t i = 0;
		for (; i < requiredColumns.size(); ++i)
		{
			if (headerSplit.size() < i || strcmp(headerSplit[i].c_str(), requiredColumns[i].c_str()) != 0)
			{
				invalidHeader = true;
				break;
			}
			setColumn(headerSplit[i]);
		}
		if (invalidHeader || (headerSplit.size() >= 9 && strcmp(headerSplit[8].c_str(), "FORMAT") != 0)) // check for the FORMAT column if it exists
		{
			throw std::runtime_error("Invalid VCF Header: " + headerString);
		}
		if (headerSplit.size() == 8)
		{
			setColumn("FORMAT");
		}
		else
		{
			setColumn(headerSplit[i]);
		}
		++i; // increament i to skip FORMAT if there is one
		for (; i < headerSplit.size(); ++i)
		{
			setColumn(headerSplit[i]);
			m_sample_names_by_column_order.emplace_back(headerSplit[i]);
			m_sample_names.emplace(headerSplit[i]);
		}
	}

	void VCFHeader::addFormatToHeader()
	{
	    m_lines.emplace_back("##FORMAT=<ID=DP_NFP,Number=1,Type=Integer,Description=\"Read count at 95 percent Smith Waterman score or above\">");
		m_lines.emplace_back("##FORMAT=<ID=DP_NP,Number=1,Type=Integer,Description=\"Read count between 90 and 94 percent Smith Waterman score\">");
		m_lines.emplace_back("##FORMAT=<ID=DP_EP,Number=1,Type=Integer,Description=\"Read count between 80 and 89 percent Smith Waterman score\">");
		m_lines.emplace_back("##FORMAT=<ID=DP_SP,Number=1,Type=Integer,Description=\"Read count between 70 and 79 percent Smith Waterman score\">");
		m_lines.emplace_back("##FORMAT=<ID=DP_LP,Number=1,Type=Integer,Description=\"Read count at 69 percent or less Smith Waterman score\">");
		m_lines.emplace_back("##FORMAT=<ID=DP_AP,Number=1,Type=Integer,Description=\"Read count for mappings which map equally well into (or out of) reference and variant. Not resolvable but valid mapping.\">");
		m_lines.emplace_back("##FORMAT=<ID=DP4_NFP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling at 95 percent Smith Waterman score or above.\">");
		m_lines.emplace_back("##FORMAT=<ID=DP4_NP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling between 90 and 94 percent Smith Waterman score.\">");
		m_lines.emplace_back("##FORMAT=<ID=DP4_EP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling between 80 and 89 percent Smith Waterman score.\">");
		m_lines.emplace_back("##FORMAT=<ID=DP4_SP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling between 70 and 79 percent Smith Waterman score.\">");
		m_lines.emplace_back("##FORMAT=<ID=DP4_UP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling at 69 percent or less Smith Waterman score.\">");
	}

	std::vector< std::string > VCFHeader::getColumnNames()
	{
		return this->m_columns;
	}

	int32_t VCFHeader::getColumnPosition(const std::string& columnTitle)
	{
		int32_t idx = -1;
		for (auto i = 0; i < m_columns.size(); ++i)
		{
			if (strcmp(columnTitle.c_str(), m_columns[i].c_str()) == 0)
			{
				idx = i;
				break;
			}
		}
		return idx;
	}

	void VCFHeader::setColumn(const std::string& column)
	{
		bool exists = false;
		for (auto i = 0; i < m_columns.size(); ++i)
		{
			if (strcmp(column.c_str(), m_columns[i].c_str()) == 0)
			{
				exists = true;
				break;
			}
		}
		if (!exists)
		{
			m_columns.emplace_back(column);
		}
	}

	std::string VCFHeader::getHeader()
	{
		std::string headerString = "";
		for (auto headerLine : this->m_lines)
		{
			headerString += headerLine + "\n";
		}
		headerString += getColumnsString() + "\n";
		return headerString;
	}

	std::string VCFHeader::getColumnsString()
	{
		std::string columnsString = "";
		for (auto i = 0; i < this->m_columns.size(); ++i)
		{
			std::string suffix = (i < this->m_columns.size() - 1) ? "\t" : "";
			columnsString += this->m_columns[i] + suffix;
		}
		return columnsString;
	}

	void VCFHeader::registerReferencePath(const std::string& referencePath)
	{
		m_reference_path = referencePath;
	}

	void VCFHeader::registerActiveSample(SampleManager::SharedPtr sampleManagerPtr)
	{
		for (auto samplePtr : sampleManagerPtr->getSamplePtrs())
		{
			if (!isActiveSampleColumnName(samplePtr->getName()))
			{
				m_sample_names_by_column_order.emplace_back(samplePtr->getName());
				this->m_sample_names.emplace(samplePtr->getName());
				this->m_active_sample_names.emplace(samplePtr->getName());
				bool isInColumns = false;
				for (auto& col : m_columns) { if (col.compare(samplePtr->getName()) == 0) { isInColumns = true; break; } }
				if (!isInColumns) { setColumn(samplePtr->getName()); }
			}
		}
	}

	bool VCFHeader::isActiveSampleColumnName(const std::string& sampleName)
	{
		return this->m_active_sample_names.find(sampleName) != this->m_active_sample_names.end();
	}

	bool VCFHeader::isSampleColumnName(const std::string& sampleName)
	{
		return this->m_sample_names.find(sampleName) != this->m_sample_names.end();
	}

	std::vector< std::string > VCFHeader::getSampleNames()
	{
		return m_sample_names_by_column_order;
	}
}
