#ifndef GWIZ_IALLELE_H
#define GWIZ_IALLELE_H

#include <core/sequence/Sequence.h>

#include <boost/noncopyable.hpp>

#include <atomic>

namespace gwiz
{
	class IAllele : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAllele > SharedPtr;
	    IAllele(Sequence::SharedPtr sequencePtr) :
		    m_sequence_ptr(sequencePtr),
			m_padding_prefix(0),
			m_padding_suffix(0),
			m_forward_count(0),
			m_reverse_count(0),
			m_total_count(0)
		{ }

		virtual ~IAllele() {}

		Sequence::SharedPtr getSequencePtr() { return this->m_sequence_ptr; }
		const char* getSequence() { return this->m_sequence_ptr->getSequence(); }
		std::string getSequenceString() { return this->m_sequence_ptr->getSequenceString(); }

		inline uint32_t getPaddingPrefix() { return this->m_padding_prefix; }
		inline uint32_t getPaddingSuffix() { return this->m_padding_suffix; }

		inline uint32_t getForwardCount() { return this->m_forward_count.load(); }
		inline uint32_t getReverseCount() { return this->m_reverse_count.load(); }
		inline uint32_t getTotalCount() { return this->m_total_count.load(); }

		inline void incrementForwardCount() { ++this->m_forward_count; ++this->m_total_count; }
		inline void incrementReverseCount() { ++this->m_reverse_count; ++this->m_total_count; }

	protected:
		Sequence::SharedPtr m_sequence_ptr;
		uint32_t m_padding_prefix;
		uint32_t m_padding_suffix;

		std::atomic< uint32_t > m_forward_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_reverse_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_total_count; // since this needs to be accessed by several threads make it atomic

	private:

	};
}

#endif //GWIZ_IALLELE_H
