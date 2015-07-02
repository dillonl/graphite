#ifndef GWIZ_EQUIVALENTALLELE_H
#define GWIZ_EQUIVALENTALLELE_H

#include "Allele.h"

namespace gwiz
{
	class EquivalentAllele : public Allele
	{
	public:
		EquivalentAllele();
		~EquivalentAllele();

		std::shared_ptr< Sequence > getSequencePtr() override { return this->m_sequence_ptr; }
		const char* getSequence() override { return this->m_sequence_ptr->getSequence(); }
		std::string getSequenceString() override { return this->m_sequence_ptr->getSequenceString(); }

		void setSequence(Sequence::SharedPtr sequence) { this->m_sequence_ptr = sequence; }
		void addAllele(IAllele::SharedPtr allelePtr) { this->m_allele_ptrs.emplace_back(allelePtr); }
	private:
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		Sequence::SharedPtr m_sequence_ptr;
	};
}

#endif //GWIZ_EQUIVALENTALLELE_H
