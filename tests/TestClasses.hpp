#ifndef GRAPHITE_TEST_TESTCLASSES_HPP
#define GRAPHITE_TEST_TESTCLASSES_HPP

class VariantTest : public graphite::Variant
{
public:
	void setAlleleCounts(std::vector< std::string > alleles, std::vector< std::tuple< uint32_t, uint32_t > >& alleleCounts)
	{
		this->m_total_allele_count = 0;
		this->m_alt_allele_ptrs.clear();
		this->m_allele_count.clear();
		if (alleles.size() != alleleCounts.size()) { throw "alleles must be matched by allele counts"; }
		for (uint32_t i = 0; i < alleles.size(); ++i)
		{
			auto allelePtr = std::make_shared< graphite::Allele >(alleles[i]);
			this->m_total_allele_count += std::get< 0 >(alleleCounts[i]) + std::get< 1 >(alleleCounts[i]);
			this->m_allele_count[alleles[i]] = alleleCounts[i];
			if (i == 0)	{ m_ref_allele_ptr = allelePtr; }
			else { m_alt_allele_ptrs.emplace_back(allelePtr); }
		}
		this->m_all_allele_ptrs.reserve(this->m_alt_allele_ptrs.size() + 1);
		this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
		this->m_all_allele_ptrs.insert(this->m_all_allele_ptrs.end(), this->m_alt_allele_ptrs.begin(), this->m_alt_allele_ptrs.end());
	}

protected:
	std::unordered_map< std::string, std::tuple< uint32_t, uint32_t > > m_allele_count;
	uint32_t m_total_allele_count;
};



#endif //GRAPHITE_TEST_TESTCLASSES_HPP
