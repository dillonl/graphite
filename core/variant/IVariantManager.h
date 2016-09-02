#ifndef GRAPHITE_IVARIANTMANAGER_H
#define GRAPHITE_IVARIANTMANAGER_H

#include "IVariantList.h"
#include "core/region/Region.h"
#include "core/util/Noncopyable.hpp"

#include <memory>

namespace graphite
{
	class IVariantManager : Noncopyable
	{
	public:
		typedef std::shared_ptr< IVariantManager > SharedPtr;
		IVariantManager() {}
		virtual ~IVariantManager() {}

		virtual IVariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr) = 0;
		virtual IVariantList::SharedPtr getCompleteVariantList() = 0;
		virtual void releaseResources() = 0;
	};
}

#endif //GRAPHITE_IVARIANTMANAGER_H
