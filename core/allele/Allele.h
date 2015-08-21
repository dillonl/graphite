#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core/sequence/SequenceManager.h"
#include "IAllele.h"
#include "AlleleMetaData.h"

namespace graphite
{
	class VCFFileReader;
	class Allele : public IAllele
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
	    Allele(Sequence::SharedPtr sequencePtr) :
		    m_forward_count(0),
			m_reverse_count(0),
		    m_sequence_ptr(sequencePtr),
			m_allele_meta_data_ptr(std::make_shared< AlleleMetaData >(0, 0))
		{
		}

	    Allele(const std::string& seqString) :
		    m_forward_count(0),
			m_reverse_count(0),
			m_sequence_ptr(SequenceManager::Instance()->getSequence(seqString)),
			m_allele_meta_data_ptr(std::make_shared< AlleleMetaData >(0, 0))
		{
		}

	    Allele(const std::string& seqString, AlleleMetaData::SharedPtr alleleMetaDataPtr) :
		    m_forward_count(0),
			m_reverse_count(0),
			m_sequence_ptr(SequenceManager::Instance()->getSequence(seqString)),
			m_allele_meta_data_ptr(alleleMetaDataPtr)
		{
		}

		~Allele()
		{
		}


		size_t getLength() override { return this->m_sequence_ptr->getLength(); }
		std::shared_ptr< Sequence > getSequencePtr() override { return this->m_sequence_ptr; }
		const char* getSequence() override { return this->m_sequence_ptr->getSequence(); }
		void setSequence(std::shared_ptr< Sequence > sequencePtr) override { this->m_sequence_ptr = sequencePtr; }
		std::string getSequenceString() override { return this->m_sequence_ptr->getSequenceString(); }
		virtual void setAlleleMetaData(AlleleMetaData::SharedPtr alleleMetaDataPtr)  override { this->m_allele_meta_data_ptr = alleleMetaDataPtr; }
		virtual AlleleMetaData::SharedPtr getAlleleMetaData() override { return this->m_allele_meta_data_ptr; }

		virtual inline uint32_t getForwardCount() override { return this->m_forward_count.load(); }
		virtual inline uint32_t getReverseCount() override { return this->m_reverse_count.load(); }
		virtual inline uint32_t getTotalCount() override { return this->getForwardCount() + this->getReverseCount(); }

		virtual inline void incrementForwardCount() override { ++this->m_forward_count; }
		virtual inline void incrementReverseCount() override { ++this->m_reverse_count; }

	protected:
		Allele() {}

		std::atomic< uint32_t > m_forward_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_reverse_count; // since this needs to be accessed by several threads make it atomic

		Sequence::SharedPtr m_sequence_ptr;
		AlleleMetaData::SharedPtr m_allele_meta_data_ptr;

	};
}

#endif //GRAPHITE_ALLELE_H
