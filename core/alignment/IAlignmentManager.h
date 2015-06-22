#ifndef GWIZ_IALIGNMENTMANAGER_H
#define GWIZ_IALIGNMENTMANAGER_H

#include "IAlignmentList.h"
#include "core/region/Region.h"

#include <boost/noncopyable.hpp>

#include <memory>

namespace gwiz
{
	class IAlignmentManager : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentManager > SharedPtr;
		IAlignmentManager() {}
		~IAlignmentManager() {}

        virtual IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr) = 0;
		virtual void releaseResources() = 0;
	};
}

#endif //GWIZ_IALIGNMENTMANAGER_H
