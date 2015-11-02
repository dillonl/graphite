#include "Variant.h"
#include "core/alignment/IAlignment.h"

#include <unordered_set>

namespace graphite
{
	Variant::Variant(position pos, const std::string& chrom, const std::string& id, const std::string& quality, const std::string& filter, IAllele::SharedPtr refAllelePtr, std::vector< IAllele::SharedPtr > altAllelePtrs) : m_position(pos), m_chrom(chrom), m_id(id), m_qual(quality), m_filter(filter), m_unmapped_to_mapped_count(0), m_mapped_to_unmapped_count(0), m_repositioned_count(0)
	{
		this->m_ref_allele_ptr = refAllelePtr;
		this->m_alt_allele_ptrs = altAllelePtrs;

		this->m_all_allele_ptrs.reserve(this->m_alt_allele_ptrs.size() + 1);
		this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
		this->m_all_allele_ptrs.insert(this->m_all_allele_ptrs.end(), this->m_alt_allele_ptrs.begin(), this->m_alt_allele_ptrs.end());
	}

	Variant::Variant() :
		m_unmapped_to_mapped_count(0), m_mapped_to_unmapped_count(0), m_repositioned_count(0)
	{
	}

	Variant::~Variant()
	{
	}

	/*
	uint32_t Variant::getTotalAlleleCount()
	{
		uint32_t totalCount = 0;
		totalCount += this->m_ref_allele_ptr->getForwardCount();
		totalCount += this->m_ref_allele_ptr->getReverseCount();
		for (auto allelePtr : getAltAllelePtrs())
		{
			totalCount += allelePtr->getForwardCount();
			totalCount += allelePtr->getReverseCount();
		}
		return totalCount;
	}

	std::string Variant::getAlleleCountString()
	{
		std::string alleleCountString = "";
		alleleCountString += std::to_string(this->m_ref_allele_ptr->getForwardCount()) + ",";
		alleleCountString += std::to_string(this->m_ref_allele_ptr->getReverseCount());
		for (auto altAllelePtr : getAltAllelePtrs())
		{
			alleleCountString += ",";
			alleleCountString += std::to_string(altAllelePtr->getForwardCount()) + ",";
			alleleCountString += std::to_string(altAllelePtr->getReverseCount());
		}
		return alleleCountString;
	}
	*/

	std::string Variant::alleleString()
	{
		std::stringstream ss;
		const auto lastIter = this->m_alt_allele_ptrs.end() - 1;
		for (auto iter = this->m_alt_allele_ptrs.begin(); iter != this->m_alt_allele_ptrs.end(); ++iter)
		{
			ss << (*iter)->getSequence();
			if (iter != lastIter) { ss << ","; } // don't add a comma on the end of the alt section
		}
		return ss.str();
	}

	uint32_t Variant::getAllelePrefixOverlapMaxCount(IAllele::SharedPtr allelePtr)
	{
		auto iter = this->m_allele_prefix_max_overlap_map.find(allelePtr);
		return (iter != this->m_allele_prefix_max_overlap_map.end()) ? iter->second : 0;
	}

	uint32_t Variant::getAlleleSuffixOverlapMaxCount(IAllele::SharedPtr allelePtr)
	{
		auto iter = this->m_allele_suffix_max_overlap_map.find(allelePtr);
		return (iter != this->m_allele_suffix_max_overlap_map.end()) ? iter->second : 0;
	}

	void Variant::setAlleleOverlapMaxCountIfGreaterThan(IAllele::SharedPtr allelePtr, std::unordered_map< IAllele::SharedPtr, uint32_t >& alleleOverlapCountMap, uint32_t overlapCount)
	{
		auto iter = alleleOverlapCountMap.find(allelePtr);
		if (iter == alleleOverlapCountMap.end() || iter->second < overlapCount)
		{
			alleleOverlapCountMap[allelePtr] = overlapCount;
		}
	}

	/*
	 * This is where allele's variant weakptr is set.
	 * Also, it is with respect to this variant that
	 * the allele prefix and suffix match is calculated.
	 */
	void Variant::processOverlappingAlleles()
	{
		std::weak_ptr< IVariant > weakPtr = shared_from_this();
		uint32_t alleleCount = 1;
		m_max_prefix_match_length = 0;
		m_max_suffix_match_length = 0;
		for (auto& allelePtr1 : this->m_all_allele_ptrs)
		{
			allelePtr1->setVariantWPtr(weakPtr);
			for (uint32_t i = alleleCount; i < this->m_all_allele_ptrs.size(); ++i)
			{
				auto allelePtr2 = this->m_all_allele_ptrs[i];
				auto tmpMaxPrefixCount = allelePtr1->getCommonPrefixSize(allelePtr2);
				auto tmpMaxSuffixCount = allelePtr1->getCommonSuffixSize(allelePtr2);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr1, this->m_allele_prefix_max_overlap_map, tmpMaxPrefixCount);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr2, this->m_allele_prefix_max_overlap_map, tmpMaxPrefixCount);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr1, this->m_allele_suffix_max_overlap_map, tmpMaxSuffixCount);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr2, this->m_allele_suffix_max_overlap_map, tmpMaxSuffixCount);
			}
			++alleleCount;
		}
	}

	void Variant::setRefAndAltAlleles(const std::string& ref, const std::vector< std::string >& alts)
	{
		this->m_all_allele_ptrs.clear();
		this->m_ref_allele_ptr = std::make_shared< Allele >(SequenceManager::Instance()->getSequence(ref.c_str()));
		this->m_all_allele_ptrs.reserve(alts.size() + 1);
		this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
		this->m_alt_allele_ptrs.clear();
		for (const auto& alt : alts)
		{
			auto sequencePtr = SequenceManager::Instance()->getSequence(alt.c_str());
			auto altAllelePtr = std::make_shared< Allele >(sequencePtr);
			this->m_alt_allele_ptrs.emplace_back(altAllelePtr);
			this->m_all_allele_ptrs.emplace_back(altAllelePtr);
		}
	}

	void Variant::incrementUnmappedToMappedCount()
	{
		++this->m_unmapped_to_mapped_count;
	}

	void Variant::incrementMappedToUnmappedCount()
	{
		++this->m_mapped_to_unmapped_count;
	}

	void Variant::incrementRepositionedCount()
	{
		++this->m_repositioned_count;
	}

	std::string Variant::getGenotype()
	{
		return "";
	}

	std::string Variant::getInfoFieldsString()
	{
		std::string infoFields = "";
		for (auto infoField : this->m_info_fields)
		{
			std::string prefix = (infoFields.size() > 0) ? "\t" : "";
			infoFields += prefix + infoField.first + "=" + infoField.second;
		}
		return (this->m_info_fields.size() > 0) ? infoFields : ".";
	}

	void Variant::printVariant(std::ostream& out, std::vector< std::shared_ptr< Sample > > samplePtrs)
	{
		out << this->m_chrom << "\t" << getPosition() << "\t.\t" << this->m_ref_allele_ptr->getSequence() << "\t" << alleleString() << "\t0\t.\t" << getInfoFieldsString() << getSampleCounts(samplePtrs) << std::endl;
	}

	std::string Variant::getSampleCounts(std::vector< Sample::SharedPtr > samplePtrs)
	{
		std::string alleleCountString = "";
		std::unordered_map< std::string, std::vector< VariantSampleContainer > > sampleNameToCounts; // counts are total, ( forward, reverse )

		std::unordered_set< std::string > sampleNameSet;
		for (auto samplePtr : samplePtrs)
		{

			if (sampleNameSet.find(samplePtr->getName()) == sampleNameSet.end())
			{
				sampleNameSet.emplace(samplePtr->getName());

				uint32_t totalCount = 0;
				std::string sampleString = "";
				for (size_t i = 0; i < m_all_allele_ptrs.size(); ++i)
				{
					auto allelePtr = m_all_allele_ptrs[i];
					uint32_t forwardCount = allelePtr->getForwardCount(samplePtr);
					uint32_t reverseCount = allelePtr->getReverseCount(samplePtr);
					std::string prefix = (i == 0) ? "" : ",";
					sampleString += prefix + std::to_string(forwardCount) + "," + std::to_string(reverseCount);
					totalCount += (forwardCount + reverseCount);
				}
				alleleCountString += "\tDP=" + std::to_string(totalCount) + ";DP4=" + sampleString;
			}
		}
		return alleleCountString;
	}

}
