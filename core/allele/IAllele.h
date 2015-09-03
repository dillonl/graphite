#ifndef GRAPHITE_IALLELE_H
#define GRAPHITE_IALLELE_H

#include "AlleleMetaData.h"

#include <boost/noncopyable.hpp>

#include <atomic>
#include <unordered_map>

namespace graphite
{
	class VCFReader;
	class Sequence;
	class IVariant;
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

		virtual uint32_t getForwardCount() = 0;
		virtual uint32_t getReverseCount() = 0;
		virtual uint32_t getTotalCount() = 0;

		virtual void incrementForwardCount() = 0;
		virtual void incrementReverseCount() = 0;

		void setVariantWPtr(std::weak_ptr< IVariant > variantWPtr) { m_variant_wptr = variantWPtr; }
		std::weak_ptr< IVariant > getVariantWPtr() { return m_variant_wptr; }

		inline SharedPtr getSharedPtr()
		{
			return shared_from_this();
		}

		virtual uint32_t getCommonPrefixSize(IAllele::SharedPtr allelePtr) = 0;
		virtual uint32_t getCommonSuffixSize(IAllele::SharedPtr allelePtr) = 0;

	protected:

			std::weak_ptr< IVariant > m_variant_wptr;

	};
}

#endif //GRAPHITE_IALLELE_H
