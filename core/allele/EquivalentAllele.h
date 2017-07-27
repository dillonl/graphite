#ifndef GRAPHITE_EQUIVALENTALLELE_H
#define GRAPHITE_EQUIVALENTALLELE_H

#include "Allele.h"

namespace graphite
{
	/*
	 * The purpose of this allele class is to contain alleles
	 * that have the same sequence but different padding and
	 * therefore can be part of multiple VCF files. The reason
	 * that the alleles are contained all within this allele
	 * is so that when it is adjudicated the adjudication can
	 * be propagated to all alleles in m_allele_ptrs.
	 */
	class EquivalentAllele : public Allele
	{
	public:
	    EquivalentAllele(const std::string& sequence) : Allele(sequence)
		{
			int x = 0;
			if (x > 0)
			{
				x = 1;
			}
		}
	    EquivalentAllele(const std::string& seqString, std::vector< IAllele::SharedPtr > allelePtrs) :
		    Allele(seqString),
			m_allele_ptrs(allelePtrs)
		{
		}
		~EquivalentAllele() {}

		void addAllele(IAllele::SharedPtr allelePtr) { this->m_allele_ptrs.emplace_back(allelePtr); }
		std::vector< IAllele::SharedPtr > getAllAlleles() { return this->m_allele_ptrs; }

		inline uint32_t getForwardCount(const std::string& sampleName, AlleleCountType alleleCountType) override
		{
			uint32_t count = 0;
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				count += allelePtr->getForwardCount(sampleName, alleleCountType);
			}
			return count;
		}
		inline uint32_t getReverseCount(const std::string& sampleName, AlleleCountType alleleCountType) override
		{
			uint32_t count = 0;
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				count += allelePtr->getReverseCount(sampleName, alleleCountType);
			}
			return count;
		}
		inline uint32_t getTotalCount(AlleleCountType alleleCountType) override
		{
			uint32_t count = 0;
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				count += allelePtr->getTotalCount(alleleCountType);
			}
			return count;
		}

		inline void incrementForwardCount(std::shared_ptr< Sample > samplePtr, AlleleCountType alleleCountType) override
		{
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				allelePtr->incrementForwardCount(samplePtr, alleleCountType);
			}
		}
		inline void incrementReverseCount(std::shared_ptr< Sample > samplePtr, AlleleCountType alleleCountType) override
		{
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				allelePtr->incrementReverseCount(samplePtr, alleleCountType);
			}
		}
		inline void incrementCount(bool isReverseStrand, std::shared_ptr< Sample > samplePtr, AlleleCountType alleleCountType) override
		{
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				allelePtr->incrementCount(isReverseStrand, samplePtr, alleleCountType);
			}
		}

		/*
		void addCommonPrefixInformation(uint32_t prefixSize, IAllele::SharedPtr allelePtrs) override
		{
		}
		void addCommonSuffixInformation(uint32_t suffixSize, IAllele::SharedPtr allelePtrs) override
		{
		}
		*/
	private:
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
	};
}

#endif //GRAPHITE_EQUIVALENTALLELE_H
