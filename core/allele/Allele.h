#ifndef GWIZ_ALLELE_H
#define GWIZ_ALLELE_H

#include "core/sequence/Sequence.h"
#include "IAllele.h"
#include "AlleleMetaData.h"

namespace gwiz
{
	class VCFFileReader;
	class Allele : public IAllele
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
	    Allele(std::shared_ptr< Sequence > sequencePtr) :
		    m_sequence_ptr(sequencePtr)
		{
		}

		Allele() = delete;
		~Allele() {}

		std::shared_ptr< Sequence > getSequencePtr() override { return this->m_sequence_ptr; }
		const char* getSequence() override { return this->m_sequence_ptr->getSequence(); }
		std::string getSequenceString() override { return this->m_sequence_ptr->getSequenceString(); }

	protected:
		std::shared_ptr< Sequence > m_sequence_ptr;

	};
}

#endif //GWIZ_ALLELE_H
