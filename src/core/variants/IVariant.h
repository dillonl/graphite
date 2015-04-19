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
            IVariant() {}
            virtual ~IVariant() {}

            virtual size_t getSmallestAlleleSize() = 0;
            virtual size_t getLargestAlleleSize() = 0;

            virtual std::string toString() = 0;
        private:
    };
}

#endif// GWIZ_IVARIANT_H
