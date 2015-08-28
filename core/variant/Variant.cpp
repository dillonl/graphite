#include "Variant.h"
#include "core/alignment/IAlignment.h"

namespace graphite
{
	Variant::Variant(position pos, const std::string& chrom, const std::string& id, const std::string& quality, const std::string& filter, IAllele::SharedPtr refAllelePtr, std::vector< IAllele::SharedPtr > altAllelePtrs) : m_position(pos), m_chrom(chrom), m_id(id), m_qual(quality), m_filter(filter)
	{
		this->m_ref_allele_ptr = refAllelePtr;
		this->m_alt_allele_ptrs = altAllelePtrs;

		this->m_all_allele_ptrs.reserve(this->m_alt_allele_ptrs.size() + 1);
		this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
		this->m_all_allele_ptrs.insert(this->m_all_allele_ptrs.end(), this->m_alt_allele_ptrs.begin(), this->m_alt_allele_ptrs.end());
	}

	Variant::Variant()
	{
	}

	Variant::~Variant()
	{
	}

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

	/*
	void Variant::processOverlappingAlleles()
	{
		std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > > sharedPrefixes;
		std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > > sharedSuffixes;
		std::unordered_map< IAllele::SharedPtr, std::unordered_map< IAllele::SharedPtr, bool > > alleleComputedChecker;
		for (auto& allelePtr1 : this->m_all_allele_ptrs)
		{
			for (auto& allelePtr2 : this->m_all_allele_ptrs)
			{
				auto iter = alleleComputedChecker.find(allelePtr2);
				if (allelePtr1.get() == allelePtr2.get() || (iter != alleleComputedChecker.end() && iter->second.find(allelePtr1) != iter->second.end())) { continue; } // don't compare identical alleles or alleles where the prefix/suffix have already been computed
				auto commonPrefixSize = allelePtr1->getCommonPrefixSize(allelePtr2);
				auto commonSuffixSize = allelePtr1->getCommonSuffixSize(allelePtr2);
				std::cout << "cps "<< commonPrefixSize << std::endl;
				std::cout << "css "<< commonSuffixSize << std::endl;
				if (commonPrefixSize > 0)
				{
					allelePtr1->addCommonPrefixInformation(commonPrefixSize, allelePtr2);
					allelePtr2->addCommonPrefixInformation(commonPrefixSize, allelePtr1);
				}
				if (commonSuffixSize > 0)
				{
					allelePtr1->addCommonSuffixInformation(commonSuffixSize, allelePtr2);
					allelePtr2->addCommonSuffixInformation(commonSuffixSize, allelePtr1);
				}
				alleleComputedChecker[allelePtr2].emplace(allelePtr1, true);
			}
		}
	}
	*/

	std::string Variant::getGenotype()
	{
		/*
		std::lock_guard< std::mutex > lock(this->m_allele_count_mutex);
		float thresholdPercent = 0.3; // if there isn't at least this percent represented in the allele count then we don't consider it
		std::vector< std::tuple< std::string, uint32_t > > potentialGenotypes;

		// get the reference count
		auto rc = this->m_allele_count[this->getRef()];
		uint32_t rCount = std::get< 0 >(rc) + std::get< 1 >(rc);
		if (rCount > 3 && rCount >= (this->m_total_allele_count * thresholdPercent)) // if the reference count meets our criteria then add it to our potential genotypes
		{
			potentialGenotypes.emplace_back(std::make_tuple("0", rCount));
		}

		for (uint32_t i = 0; i < this->m_alt.size(); ++i)
		{
			auto ac = this->m_allele_count[this->m_alt[i]];
			uint32_t aCount = std::get< 0 >(ac) + std::get< 1 >(ac);
			if (aCount > 3 && aCount >= (this->m_total_allele_count * thresholdPercent)) // if the alt count meets our criteria then add it to our potential genotypes
			{
				uint32_t alleleIndex = i + 1; // we add 1 because 0 is the ref and m_alt is 0 based
				potentialGenotypes.emplace_back(std::make_tuple(std::to_string(alleleIndex), aCount));
			}
		}

		// sort the genotypes in decending order. Then we just worry about the first 2 indices (if they exist)
		std::sort(potentialGenotypes.begin(), potentialGenotypes.end(), [](const std::tuple< std::string, uint32_t >& a, const std::tuple< std::string, uint32_t >& b)
				  {
					  return std::get< 1 >(a) > std::get< 1 >(b);
				  });
		std::string genotype;
		// based on the size of the genotypes we generate the genotype
		switch (potentialGenotypes.size())
		{
		case 0: // there isn't a decernible genotype
			genotype = "./.";
			break;
		case 1: // if there is only one then we know this is a homozygous
			genotype = std::get< 0 >(potentialGenotypes[0]) + "/" + std::get< 0 >(potentialGenotypes[0]);
			break;
		default: // at least 2 elements we know (since it's sorted) that the 2 highest must be at index 0 and 1 respectively (heterozygous)
			genotype = std::get< 0 >(potentialGenotypes[0]) + "/" + std::get< 0 >(potentialGenotypes[1]);
		}
		return genotype;
		*/
		return "";
	}

	void Variant::printVariant(std::ostream& out)
	{
		out << this->m_chrom << "\t" << getPosition() << "\t.\t" << this->m_ref_allele_ptr->getSequence() << "\t" << alleleString() << "\t0\t.\tDP=" << getTotalAlleleCount() << ";DP4=" << getAlleleCountString() << std::endl;
	}

}
