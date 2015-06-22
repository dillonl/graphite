#ifndef GWIZ_IVARIANTMANAGER_H
#define GWIZ_IVARIANTMANAGER_H

#include "IVariantList.h"
#include "core/region/Region.h"

#include <boost/noncopyable.hpp>

#include <memory>

namespace gwiz
{
	class IVariantManager : boost::noncopyable
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

#endif //GWIZ_IVARIANTMANAGER_H
