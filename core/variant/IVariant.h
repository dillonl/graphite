#ifndef GRAPHITE_IVARIANT_H
#define GRAPHITE_IVARIANT_H

#include "core/util/Types.h"
#include "core/allele/IAllele.h"
#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"
#include "IHeader.h"

#include <vector>
#include <memory>
#include <unordered_set>

namespace graphite
{
	class IAllele;
	class IAlignment;
	class Sample;
    class IVariant : private Noncopyable, public std::enable_shared_from_this< IVariant >
    {
        public:
		    typedef std::shared_ptr<IVariant> SharedPtr;
			typedef std::shared_ptr<IVariant> WeakPtr;
            IVariant()
            {
            }
            virtual ~IVariant() {}

			virtual position getPosition() = 0;
			virtual std::string getChrom() const = 0;
			virtual IAllele::SharedPtr getRefAllelePtr() = 0;
			virtual std::vector< IAllele::SharedPtr > getAltAllelePtrs() = 0;
			virtual void processOverlappingAlleles() = 0; // set all allele variantwptrs to be this
			virtual uint32_t getAllelePrefixOverlapMaxCount(IAllele::SharedPtr allelePtr) = 0;
			virtual uint32_t getAlleleSuffixOverlapMaxCount(IAllele::SharedPtr allelePtr) = 0;
			virtual std::string getVariantLine(IHeader::SharedPtr headerPtr) = 0;
			virtual bool shouldSkip() = 0;
			virtual void setSkip(bool) = 0;
			virtual std::vector< Region::SharedPtr > getRegions() = 0;
			virtual bool doesOverlap(IVariant::SharedPtr variantPtr) = 0;
			virtual uint32_t getReferenceSize() = 0;
			virtual void addRegion(Region::SharedPtr regionPtr) = 0;
    };
}

#endif// GRAPHITE_IVARIANT_H
