#ifndef GWIZ_IVARIANT_H
#define GWIZ_IVARIANT_H

#include <boost/noncopyable.hpp>

#include <memory>

namespace gwiz
{
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

            virtual size_t getSmallestAlleleSize() = 0;
            virtual size_t getLargestAlleleSize() = 0;
            uint32_t getVariantID() { return this->m_variant_id; }

			virtual void printVariant(std::ostream& out) = 0;
        private:
            uint32_t m_variant_id;
    };
}

#endif// GWIZ_IVARIANT_H
