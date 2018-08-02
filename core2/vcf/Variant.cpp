#include "Variant.h"
#include "core2/util/Utility.h"
#include "core2/util/Types.h"

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
		std::string vcfLine;

		for (uint32_t i = 0; i < STANDARD_VCF_COLUMN_NAMES.size(); ++i)
		{
			if (i > 0)
			{
				vcfLine += "\t";
			}
			std::string formatStr = "";
			if (STANDARD_VCF_COLUMN_NAMES[i].compare("FORMAT") == 0)
			{
				formatStr = ":DP_NFP:DP4_NFP:DP_NP:DP4_NP:DP_EP:DP4_EP:DP_SP:DP4_SP:DP_LP:DP4_LP:DP_AP:DP4_AP";
			}
			vcfLine += m_columns[STANDARD_VCF_COLUMN_NAMES[i]] + formatStr;
		}
		auto samplePtrs = this->m_vcf_writer_ptr->getSamplePtrs();
		for (uint32_t i = 0; i < samplePtrs.size(); ++i)
		{
			std::string sep = (!m_columns[samplePtrs[i]->getName()].empty()) ? ":" : "";
			vcfLine += "\t" + m_columns[samplePtrs[i]->getName()] + sep + getGraphiteCounts(samplePtrs[i]->getName());
			// vcfLine += "\t" + getGraphiteCounts(samplePtrs[i]->getName());
		}

		this->m_vcf_writer_ptr->writeLine(vcfLine);
	}

	void Variant::parseColumns()
	{
		auto samplePtrs = this->m_vcf_writer_ptr->getSamplePtrs();
		std::vector< std::string > columns;
		split(this->m_variant_line, '\t', columns);
		if (columns.size() != (STANDARD_VCF_COLUMN_NAMES.size() + samplePtrs.size()))
		{
			std::cout << "Invalid VCF at line: " << this->m_variant_line << std::endl;
			exit(EXIT_FAILURE);
		}
		for (uint32_t i = 0; i < STANDARD_VCF_COLUMN_NAMES.size(); ++i)
		{
			auto standardColumn = STANDARD_VCF_COLUMN_NAMES[i];
			m_columns[standardColumn] = columns[i];
		}
		this->m_chrom = m_columns["#CHROM"];
		this->m_position = stoi(m_columns["POS"]);
		for (uint32_t i = 0; i < samplePtrs.size(); ++i)
		{
			size_t columnIdx = STANDARD_VCF_COLUMN_NAMES.size() + i;
			m_columns[samplePtrs[i]->getName()] = columns[columnIdx];
		}
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

	/*
	std::vector< Node::SharedPtr > Variant::getReferenceNodePtrs()
	{
		std::vector< Node::SharedPtr > referenceNodePtrs;
		std::unordered_set< Node::SharedPtr > altNodePtrs;
		for (auto altAllelePtr : this->m_alternate_allele_ptrs)
		{
			altNodePtrs.emplace(altAllelePtr->getNodePtr());
		}
		Node::SharedPtr refNodePtr = this->m_alternate_allele_ptrs[0]->getNodePtr()->getReferenceInNode();
		// std::vector< Node::SharedPtr > inNodePtrs = this->m_alternate_allele_ptrs[0]->getNodePtr()->getInNodes();
		// Node::SharedPtr refNodePtr = inNodePtrs[0];
		do
		{
			for (auto nodePtr : refNodePtr->getOutNodes())
			{
				if (nodePtr->getAlleleType() == Node::ALLELE_TYPE::REF)
				{
					refNodePtr = nodePtr;
					break;
				}
			}
			for (auto nodePtr : refNodePtr->getInNodes())
			{
				if (altNodePtrs.count(nodePtr) > 0)
				{
					altNodePtrs.erase(nodePtr);
				}
			}
			referenceNodePtrs.emplace_back(refNodePtr);
		} while (!altNodePtrs.empty());
		return referenceNodePtrs;
	}
	*/

	std::string Variant::getGraphiteCounts(const std::string& sampleName)
	{
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
			for (auto allelePtr : this->m_alternate_allele_ptrs)
			{
				tmpCountsString += ",";
				forwardScoreCount = allelePtr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, true);
				reverseScoreCount = allelePtr->getScoreCountFromAlleleCountType(sampleName, alleleCountType, false);
				tmpCountsString += std::to_string(forwardScoreCount.size()) + "," + std::to_string(reverseScoreCount.size());
				/*
				auto tmpCount123 = forwardScoreCount.size() + reverseScoreCount.size();
				if (tmpCount123 > 0)
				{
					std::cout << "alt count found" << std::endl;
				}
				*/
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
