#ifndef GWIZ_IVARIANT_H
#define GWIZ_IVARIANT_H

#include "utils/NonCopyable.h"

namespace gwiz
{
    class IVariant : private noncopyable
    {
        public:
            IVariant() {}
            virtual ~IVariant() {}
        private:
            IVariant(const IVariant& other) = delete;
            IVariant& operator=(const IVariant& other) = delete;
    };
}

#endif// GWIZ_IVARIANT_H
