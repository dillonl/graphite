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
	    EquivalentAllele(Sequence::SharedPtr sequence) : Allele(sequence) {}
	    EquivalentAllele(const std::string& seqString, std::vector< IAllele::SharedPtr > allelePtrs) :
		    Allele(seqString),
			m_allele_ptrs(allelePtrs)
		{
		}
		~EquivalentAllele() {}

		void addAllele(IAllele::SharedPtr allelePtr) { this->m_allele_ptrs.emplace_back(allelePtr); }
		std::vector< IAllele::SharedPtr > getAllAlleles() { return this->m_allele_ptrs; }

		inline uint32_t getForwardCount() override
		{
			uint32_t count = 0;
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				count += allelePtr->getForwardCount();
			}
			return count;
		}
		inline uint32_t getReverseCount() override
		{
			uint32_t count = 0;
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				count += allelePtr->getReverseCount();
			}
			return count;
		}
		inline uint32_t getTotalCount() override
		{
			uint32_t count = 0;
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				count += allelePtr->getForwardCount();
				count += allelePtr->getReverseCount();
			}
			return count;
		}

		inline void incrementForwardCount() override
		{
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				allelePtr->incrementForwardCount();
			}
		}
		inline void incrementReverseCount() override
		{
			for (auto& allelePtr : this->m_allele_ptrs)
			{
				allelePtr->incrementReverseCount();
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
