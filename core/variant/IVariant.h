#ifndef GWIZ_IVARIANT_H
#define GWIZ_IVARIANT_H

#include "core/util/Types.h"
#include "core/allele/IAllele.h"

#include <boost/noncopyable.hpp>

#include <vector>
#include <memory>

namespace gwiz
{
	class IAllele;
	class IAlignment;
    class IVariant : private boost::noncopyable
    {
        public:
            typedef std::shared_ptr<IVariant> SharedPtr;
            IVariant()
            {
                static uint32_t idCounter = 0;
                this->m_variant_id = ++idCounter;
            }
            virtual ~IVariant() {}

            uint32_t getVariantID() { return this->m_variant_id; }

			virtual position getPosition() = 0;
			virtual std::string getChrom() const = 0;
			virtual IAllele::SharedPtr getRefAllelePtr() = 0;
			virtual std::vector< IAllele::SharedPtr > getAltAllelePtrs() = 0;
			virtual size_t getSmallestAlleleSize() = 0;
            virtual size_t getLargestAlleleSize() = 0;
			virtual void printVariant(std::ostream& out) = 0;

			// things to be removed
			virtual void addPotentialAlignment(const std::shared_ptr< IAlignment > alignmentPtr) = 0;
        private:
            uint32_t m_variant_id;
    };
}

#endif// GWIZ_IVARIANT_H
