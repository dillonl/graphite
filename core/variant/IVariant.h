#ifndef GRAPHITE_IVARIANT_H
#define GRAPHITE_IVARIANT_H

#include "core/util/Types.h"
#include "core/allele/IAllele.h"
#include "core/util/Noncopyable.hpp"

#include <vector>
#include <memory>

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
			virtual void incrementUnmappedToMappedCount() = 0;
			virtual void incrementMappedToUnmappedCount() = 0;
			virtual void incrementRepositionedCount() = 0;
			virtual void printVariant(std::ostream& out, std::vector< std::shared_ptr< Sample > > samplePtrs) = 0;
			virtual std::string getVariantLine(std::vector< std::shared_ptr< Sample > > samplePtrs) = 0;
			virtual size_t getMaxAlleleSize() = 0;
			virtual bool shouldSkip() = 0;
			virtual void setSkip(bool) = 0;
    };
}

#endif// GRAPHITE_IVARIANT_H
