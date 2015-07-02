#ifndef GWIZ_IALLELE_H
#define GWIZ_IALLELE_H

#include "AlleleMetaData.h"

#include <boost/noncopyable.hpp>

#include <atomic>

namespace gwiz
{
	class VCFReader;
	class Sequence;
	class IAllele : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAllele > SharedPtr;
	    IAllele() { }

		virtual ~IAllele() {}

		virtual IAllele::SharedPtr copyAllele() = 0;
		virtual std::shared_ptr< Sequence > getSequencePtr() = 0;
		virtual const char* getSequence() = 0;
		virtual std::string getSequenceString() = 0;

		inline uint32_t getForwardCount() { return this->m_forward_count.load(); }
		inline uint32_t getReverseCount() { return this->m_reverse_count.load(); }
		inline uint32_t getTotalCount() { return this->m_total_count.load(); }

		inline void incrementForwardCount() { ++this->m_forward_count; ++this->m_total_count; }
		inline void incrementReverseCount() { ++this->m_reverse_count; ++this->m_total_count; }

	protected:

		std::atomic< uint32_t > m_forward_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_reverse_count; // since this needs to be accessed by several threads make it atomic
		std::atomic< uint32_t > m_total_count; // since this needs to be accessed by several threads make it atomic

	private:

	};
}

#endif //GWIZ_IALLELE_H
