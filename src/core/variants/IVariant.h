#ifndef GWIZ_IVARIANT_H
#define GWIZ_IVARIANT_H

#include "utils/NonCopyable.h"

#include <memory>

namespace gwiz
{
    class IVariant : private noncopyable
    {
        public:
            typedef std::shared_ptr<IVariant> SharedPtr;
            IVariant() {}
            virtual ~IVariant() {}
        private:
            IVariant(const IVariant& other) = delete;
            IVariant& operator=(const IVariant& other) = delete;
    };
}

#endif// GWIZ_IVARIANT_H
