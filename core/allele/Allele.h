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
	    Allele(Sequence::SharedPtr sequencePtr) :
		    m_sequence_ptr(sequencePtr)
		{
		}

		Allele() = delete;
		~Allele() {}


		std::shared_ptr< Sequence > getSequencePtr() override { return this->m_sequence_ptr; }
		const char* getSequence() override { return this->m_sequence_ptr->getSequence(); }
		void setSequence(std::shared_ptr< Sequence > sequencePtr) override { this->m_sequence_ptr = sequencePtr; }
		std::string getSequenceString() override { return this->m_sequence_ptr->getSequenceString(); }
		virtual void setAlleleMetaData(AlleleMetaData::SharedPtr alleleMetaDataPtr)  override { this->m_allele_meta_data_ptr = alleleMetaDataPtr; }
		virtual AlleleMetaData::SharedPtr getAlleleMetaData() override { return this->m_allele_meta_data_ptr; }

	protected:
		Sequence::SharedPtr m_sequence_ptr;
		AlleleMetaData::SharedPtr m_allele_meta_data_ptr;
	};
}

#endif //GWIZ_ALLELE_H
