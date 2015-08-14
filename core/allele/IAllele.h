#ifndef GRAPHITE_IALLELE_H
#define GRAPHITE_IALLELE_H

#include "AlleleMetaData.h"

#include <boost/noncopyable.hpp>

#include <atomic>

namespace graphite
{
	class VCFReader;
	class Sequence;
	class IAllele : private boost::noncopyable, public std::enable_shared_from_this< IAllele >
	{
	public:
		typedef std::shared_ptr< IAllele > SharedPtr;
	    IAllele()
		{
		}

		virtual ~IAllele() {}
		virtual size_t getLength() = 0;
		virtual std::shared_ptr< Sequence > getSequencePtr() = 0;
		virtual const char* getSequence() = 0;
		virtual std::string getSequenceString() = 0;
		virtual void setSequence(std::shared_ptr< Sequence > sequencePtr) = 0;
		virtual void setAlleleMetaData(AlleleMetaData::SharedPtr alleleMetaDataPtr) = 0;
		virtual AlleleMetaData::SharedPtr getAlleleMetaData() = 0;

		virtual inline uint32_t getForwardCount() = 0;
		virtual inline uint32_t getReverseCount() = 0;
		virtual inline uint32_t getTotalCount() = 0;

		virtual inline void incrementForwardCount() = 0;
		virtual inline void incrementReverseCount() = 0;

		inline SharedPtr getSharedPtr()
		{
			return shared_from_this();
		}

	protected:

	private:

	};
}

#endif //GRAPHITE_IALLELE_H
