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

		virtual uint32_t getForwardCount() override { return this->m_forward_count.load(); }
		virtual uint32_t getReverseCount() override { return this->m_reverse_count.load(); }
		virtual uint32_t getTotalCount() override { return this->getForwardCount() + this->getReverseCount(); }

		virtual void incrementForwardCount() override { ++this->m_forward_count; }
		virtual void incrementReverseCount() override { ++this->m_reverse_count; }

		void setSharedAllelePrefixAndSuffix(std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > > sharedPrefixes, std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > > sharedSuffixes)
		{
		}

		size_t getCommonPrefixSize(IAllele::SharedPtr allelePtr) override
		{
			size_t commonIndex = 0;
			if (this->getLength() < allelePtr->getLength())
			{
				commonIndex = this->getSequenceString().find_first_not_of(allelePtr->getSequenceString());
			}
			else
			{
				commonIndex = allelePtr->getSequenceString().find_first_not_of(this->getSequenceString());
			}
			// note: if nothing is found then the commonindex is -1 (when we return we add 1 so we return 0. Excellent!!!)
			return commonIndex + 1;
		}

		size_t getCommonSuffixSize(IAllele::SharedPtr allelePtr) override
		{
			size_t commonIndex = 0;
			if (this->getLength() < allelePtr->getLength())
			{
				commonIndex = this->getSequenceString().find_last_not_of(allelePtr->getSequenceString());
			}
			else
			{
				commonIndex = allelePtr->getSequenceString().find_last_not_of(this->getSequenceString());
			}
			// note: if nothing is found then the commonindex is -1 (when we return we add 1 so we return 0. Excellent!!!)
			return commonIndex + 1;
		}

		void addCommonPrefixInformation(uint32_t prefixSize, IAllele::SharedPtr allelePtr) override
		{
			addCommonInformationHelper(m_shared_prefixes, prefixSize, allelePtr);
		}

		void addCommonSuffixInformation(uint32_t suffixSize, IAllele::SharedPtr allelePtr) override
		{
			addCommonInformationHelper(m_shared_suffixes, suffixSize, allelePtr);
		}

	protected:
		Allele() {}

		void addCommonInformationHelper(std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > >& sharedMap, uint32_t sharedSize, IAllele::SharedPtr allelePtr)
		{
			auto iter = sharedMap.find(sharedSize);
			if (iter == sharedMap.end())
			{
				sharedMap[sharedSize] = { allelePtr };
			}
			else
			{
				iter->second.emplace_back(allelePtr);
			}
		}

		std::atomic< uint32_t > m_forward_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_reverse_count; // since this needs to be accessed by several threads make it atomic

		Sequence::SharedPtr m_sequence_ptr;
		AlleleMetaData::SharedPtr m_allele_meta_data_ptr;

		std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > > m_shared_prefixes;
		std::unordered_map< uint32_t, std::vector< IAllele::SharedPtr > > m_shared_suffixes;

	};
}

#endif //GRAPHITE_ALLELE_H
