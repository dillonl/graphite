#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core/sequence/SequenceManager.h"
#include "IAllele.h"
#include "AlleleMetaData.h"
#include "core/alignment/IAlignment.h"

namespace graphite
{
	class VCFFileReader;
	class IAlignmentReader;
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

		virtual uint32_t getForwardCount(std::shared_ptr< Sample > samplePtr) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto iter = m_alignment_sample_forward_count_map.find(samplePtr);
			uint32_t totalCount = (iter == m_alignment_sample_forward_count_map.end()) ? 0 : iter->second;
			return totalCount;
		}
		virtual uint32_t getReverseCount(std::shared_ptr< Sample > samplePtr) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto iter = m_alignment_sample_reverse_count_map.find(samplePtr);
			uint32_t totalCount = (iter == m_alignment_sample_reverse_count_map.end()) ? 0 : iter->second;
			return totalCount;
		}
		virtual uint32_t getTotalCount() override
		{
			uint32_t totalCount = 0;
			for (auto iter : m_alignment_sample_forward_count_map)
			{
				totalCount += iter.second;
			}
			for (auto iter : m_alignment_sample_reverse_count_map)
			{
				totalCount += iter.second;
			}
			return totalCount;
		}

		virtual void incrementForwardCount(std::shared_ptr< IAlignment > alignmentPtr) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto samplePtr = alignmentPtr->getSample();
			auto iter = m_alignment_sample_forward_count_map.find(samplePtr);
			uint32_t count = 1;
			if (iter != m_alignment_sample_forward_count_map.end())
			{
				count = iter->second + 1;
			}
			m_alignment_sample_forward_count_map.emplace(samplePtr, count);
			/* ++this->m_forward_count; */
		}
		virtual void incrementReverseCount(std::shared_ptr< IAlignment > alignmentPtr) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto samplePtr = alignmentPtr->getSample();
			auto iter = m_alignment_sample_reverse_count_map.find(samplePtr);
			uint32_t count = 1;
			if (iter != m_alignment_sample_reverse_count_map.end())
			{
				count = iter->second + 1;
			}
			m_alignment_sample_reverse_count_map.emplace(samplePtr, count);
			/* ++this->m_reverse_count; */
		}

		uint32_t getCommonPrefixSize(IAllele::SharedPtr allelePtr) override
		{
			uint32_t commonIndex = 0;
			uint32_t minSize = (this->getLength() < allelePtr->getLength()) ? this->getLength() : allelePtr->getLength();
			for (uint32_t i = 0; i < minSize; ++i) { if (this->getSequence()[i] != allelePtr->getSequence()[i]) { return i; } }
			return minSize;
		}

		uint32_t getCommonSuffixSize(IAllele::SharedPtr allelePtr) override
		{
			uint32_t commonIndex = 0;
			uint32_t minSize = (this->getLength() < allelePtr->getLength()) ? this->getLength() : allelePtr->getLength();
			for (uint32_t i = 0; i < minSize; ++i) { if (this->getSequence()[(this->getLength() - 1) - i] != allelePtr->getSequence()[(allelePtr->getLength() - 1) - i]) { return i; } }
			return minSize;
		}

	protected:
		Allele() {}

		std::atomic< uint32_t > m_forward_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_reverse_count; // since this needs to be accessed by several threads make it atomic

		Sequence::SharedPtr m_sequence_ptr;
		AlleleMetaData::SharedPtr m_allele_meta_data_ptr;
		std::mutex m_alignment_count_mutex;
		std::unordered_map< std::shared_ptr< Sample >, uint32_t > m_alignment_sample_forward_count_map;
		std::unordered_map< std::shared_ptr< Sample >, uint32_t > m_alignment_sample_reverse_count_map;

	};
}

#endif //GRAPHITE_ALLELE_H
