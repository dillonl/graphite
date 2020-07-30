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
		if (this->m_vcf_writer_ptr->getSaveSupportingReadInfo())
		{
			std::ofstream* outPtr = this->m_vcf_writer_ptr->getSupportingReadOutputStream();
			std::string token = "\t";
			std::string variantInfo = this->m_chrom + token + std::to_string(this->m_position);
			for (auto refAlleleSupportingReadInfo : this->m_reference_allele_ptr->getSupportingReadInfoPtrs())
			{
				(*outPtr) << variantInfo << token << this->m_reference_allele_ptr->getSequence() << token << refAlleleSupportingReadInfo->toString(token) << std::endl;
			}
			for (auto altAllelePtr : this->m_alternate_allele_ptrs)
			{
				// (*outPtr) << "" << this->m_chrom << "\t" << this->m_position << "\t" << altAllelePtr->getSequence() << std::endl;
				for (auto altAlleleSupportingReadInfo : altAllelePtr->getSupportingReadInfoPtrs())
				{
					// (*outPtr) << altAlleleSupportingReadInfo->toString("\t") << std::endl;
					(*outPtr) << variantInfo << token << altAllelePtr->getSequence() << token << altAlleleSupportingReadInfo->toString(token) << std::endl;
				}
			}
		}

		this->m_vcf_writer_ptr->writeLine(vcfLine);
	}

	void Variant::processSampleColumns()
	{
		size_t colonCount = std::count(m_columns["FORMAT"].begin(), m_columns["FORMAT"].end(), ':');
		if (m_columns["FORMAT"].size() > 0)
		{
			m_columns["FORMAT"] += ":";
		}
		m_columns["FORMAT"] += "DP_NFP:DP4_NFP:DP_NP:DP4_NP:DP_EP:DP4_EP:DP_SP:DP4_SP:DP_LP:DP4_LP:DP_AP:DP2_AP:SEM";

		for (auto& sampleName : this->m_vcf_writer_ptr->getSampleNames())
		{
			std::string graphiteCounts = (this->m_vcf_writer_ptr->isSampleNameInBam(sampleName)) ? getSampleCounts(sampleName) : m_blank_graphite_format;
			std::string samplePrefix = "";
			if (!this->m_vcf_writer_ptr->isSampleNameInOriginalVCF(sampleName))
			{
				for (auto i = 0; i < colonCount; ++i)
				{
					samplePrefix += ".:";
				}
				if (samplePrefix.size() == 0 && m_columns["FORMAT"].size() > 0) // this is the scenario where there is only one entry in the format field (colonCount == 0)
				{
					samplePrefix += ".:";
				}
				m_columns[sampleName] += samplePrefix;
			}
			m_columns[sampleName] += graphiteCounts;
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
		std::string semanticString = "{";
		std::string graphiteCountsString = "";
		AlleleCountType alleleCountType = AlleleCountType::NinteyFivePercent;
		while (alleleCountType != AlleleCountType::EndEnum)
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
			uint32_t ambiguousForwardCount = 0;
			uint32_t ambiguousReverseCount = 0;
			if (alleleCountType == AlleleCountType::Ambiguous)
			{
				ambiguousForwardCount = forwardScoreCount.size();
				ambiguousReverseCount = reverseScoreCount.size();
			}
			for (auto allelePtr : this->m_alternate_allele_ptrs)
			{
				forwardScoreCount = allelePtr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, true);
				reverseScoreCount = allelePtr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, false);
				if (alleleCountType == AlleleCountType::Ambiguous)
				{
					ambiguousForwardCount += forwardScoreCount.size();
					ambiguousReverseCount += reverseScoreCount.size();
				}
				tmpCountsString += ",";
				tmpCountsString += std::to_string(forwardScoreCount.size()) + "," + std::to_string(reverseScoreCount.size());
				totalCounter += forwardScoreCount.size() + reverseScoreCount.size();
				for (auto semanticIter : allelePtr->getSemanticLocations())
				{
					if (this->m_printed_semantics.find(semanticIter.first) != this->m_printed_semantics.end())
					{
						continue;
					}
					this->m_printed_semantics.emplace(semanticIter.first);
					semanticString += std::to_string(semanticIter.first) + "(";
					int semanticCount = 0;
					for (auto semanticSequenceIter : semanticIter.second)
					{
						std::vector< std::string > alleles;
						split(semanticSequenceIter, ':', alleles);
						std::string sequences = alleles[0] + "=>" + alleles[1];
						semanticString += (semanticCount++ == 0) ? sequences : "," + sequences;
					}
					semanticString += ")";
				}
			}
			if (!graphiteCountsString.empty())
			{
				graphiteCountsString += ":";
			}
			if (alleleCountType == AlleleCountType::Ambiguous)
			{
				tmpCountsString = std::to_string(ambiguousForwardCount) + "," + std::to_string(ambiguousReverseCount);
			}
			graphiteCountsString += std::to_string(totalCounter) + ":" + tmpCountsString;
			alleleCountType = (AlleleCountType)((uint32_t)alleleCountType + 1);
		}
		if (semanticString.size() == 1) // if we didn't set any semantic alleles above
		{
			semanticString = ".";
		}
		else
		{
			semanticString += "}";
		}
		graphiteCountsString += ":" + semanticString;
		return graphiteCountsString;
	}
}
