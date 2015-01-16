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
    };
}

#endif// GWIZ_IVARIANT_H
