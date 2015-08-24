#ifndef GRAPHITE_IVARIANT_H
#define GRAPHITE_IVARIANT_H

#include "core/util/Types.h"
#include "core/allele/IAllele.h"

#include <boost/noncopyable.hpp>

#include <vector>
#include <memory>

namespace graphite
{
	class IAllele;
	class IAlignment;
    class IVariant : private boost::noncopyable
    {
        public:
            typedef std::shared_ptr<IVariant> SharedPtr;
            IVariant()
            {
            }
            virtual ~IVariant() {}

			virtual position getPosition() = 0;
			virtual std::string getChrom() const = 0;
			virtual IAllele::SharedPtr getRefAllelePtr() = 0;
			virtual std::vector< IAllele::SharedPtr > getAltAllelePtrs() = 0;
			virtual void processOverlappingAlleles() = 0;
			virtual void printVariant(std::ostream& out) = 0;
    };
}

#endif// GRAPHITE_IVARIANT_H
