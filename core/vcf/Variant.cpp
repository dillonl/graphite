#include "Variant.h"
#include "core/util/Utility.h"
#include "core/util/Types.h"

#include <algorithm>
#include <iostream>

namespace graphite
{

	Variant::Variant(const std::string& variantLine, VCFWriter::SharedPtr vcfWriterPtr) :
		m_variant_line(variantLine),
		m_vcf_writer_ptr(vcfWriterPtr),
		m_skip_adjudication(false)
	{
		parseColumns();
		setAlleles();
	}

	Variant::~Variant()
	{
	}

	void Variant::writeVariant()
	{
		std::string vcfLine = "";
		auto columnNames = m_vcf_writer_ptr->getColumnNames();
		processSampleColumns();
		for (auto i = 0; i < columnNames.size(); ++i)
		{
			if (i > 0)
			{
				vcfLine += "\t";
			}
			auto columnName = columnNames[i];
			vcfLine += m_columns[columnName];
		}
		this->m_vcf_writer_ptr->writeLine(vcfLine);
	}

	void Variant::processSampleColumns()
	{
		std::string formatSpacing = "";
		if (m_columns["FORMAT"].size() > 0)
		{
			formatSpacing = ":";
		}
		m_columns["FORMAT"] += formatSpacing + "DP_NFP:DP4_NFP:DP_NP:DP4_NP:DP_EP:DP4_EP:DP_SP:DP4_SP:DP_LP:DP4_LP:DP_AP:DP4_AP";

		for (auto& sampleName : this->m_vcf_writer_ptr->getSampleNames())
		{
			if (this->m_vcf_writer_ptr->isSampleNameInBam(sampleName))
			{
				m_columns[sampleName] += formatSpacing + getSampleCounts(sampleName);
			}
			else
			{
				m_columns[sampleName] += formatSpacing + m_blank_graphite_format;
			}
		}
	}

	void Variant::parseColumns()
	{
		m_columns.clear();
		std::vector< std::string > columns;
		split(this->m_variant_line, '\t', columns);
		auto columnNames = this->m_vcf_writer_ptr->getColumnNames();
		int columnDiff = columns.size() - columnNames.size();
		for (auto i = 0; i < columnNames.size(); ++i)
		{
			if (columnNames[i].compare("FORMAT") == 0 && this->m_vcf_writer_ptr->getBlankFormatStringPtr() == nullptr)
			{
				std::string formatCol = "";
				if (columns.size() < i)
				{
					formatCol = columns[i];
				}
				size_t n = std::count(formatCol.begin(), formatCol.end(), ':');
				std::string blankSampleFormat = (formatCol.size() > 0) ? "." : "";
				for (int i = 0; i < n; ++i) { blankSampleFormat += ":."; }
				this->m_vcf_writer_ptr->setBlankFormatString(blankSampleFormat);
			}
			if (i > columns.size() - 1)
			{
				m_columns.emplace(columnNames[i], *this->m_vcf_writer_ptr->getBlankFormatStringPtr());
			}
			else
			{
				m_columns.emplace(columnNames[i], columns[i]);
			}
		}

		this->m_chrom = m_columns["#CHROM"];
		this->m_position = stoi(m_columns["POS"]);
	}

	void Variant::setAlleles()
	{
		this->m_reference_allele_ptr = std::make_shared< Allele >(m_columns["REF"]);
		std::vector< std::string > alts;
		if (m_columns["ALT"].find(",") != std::string::npos)
		{
			split(m_columns["ALT"], ',', alts);
		}
		else
		{
			alts.emplace_back(m_columns["ALT"]);
		}
		this->m_alternate_allele_ptrs.clear();
		for (auto alt : alts)
		{
			auto altAllelePtr = std::make_shared< Allele >(alt);
			this->m_alternate_allele_ptrs.emplace_back(altAllelePtr);
		}
	}

	std::string Variant::getSampleCounts(const std::string& sampleName)
	{
		std::string graphiteCountsString = "";
		AlleleCountType alleleCountType = AlleleCountType::NinteyFivePercent;
		while (alleleCountType != AlleleCountType::EndEnum)
		// for (auto alleleCountType : AllAlleleCountTypes)
		{
			uint32_t totalCounter = 0;
			std::unordered_set< std::string > forwardScoreCount;
			std::unordered_set< std::string > reverseScoreCount;
			std::unordered_set< std::string > forwardScoreCountTmp = this->m_reference_allele_ptr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, true);
			std::unordered_set< std::string > reverseScoreCountTmp = this->m_reference_allele_ptr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, false);
			forwardScoreCount.insert(forwardScoreCountTmp.begin(), forwardScoreCountTmp.end());
			reverseScoreCount.insert(reverseScoreCountTmp.begin(), reverseScoreCountTmp.end());
			totalCounter += forwardScoreCount.size() + reverseScoreCount.size();
			std::string tmpCountsString = std::to_string(forwardScoreCount.size()) + "," + std::to_string(reverseScoreCount.size());
			for (auto allelePtr : this->m_alternate_allele_ptrs)
			{
				tmpCountsString += ",";
				forwardScoreCount = allelePtr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, true);
				reverseScoreCount = allelePtr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, false);
				tmpCountsString += std::to_string(forwardScoreCount.size()) + "," + std::to_string(reverseScoreCount.size());
				totalCounter += forwardScoreCount.size() + reverseScoreCount.size();
			}
			if (!graphiteCountsString.empty())
			{
				graphiteCountsString += ":";
			}
			graphiteCountsString += std::to_string(totalCounter) + ":" + tmpCountsString;
			alleleCountType = (AlleleCountType)((uint32_t)alleleCountType + 1);
		}

		return graphiteCountsString;
	}
}
