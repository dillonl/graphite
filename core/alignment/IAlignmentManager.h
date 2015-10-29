#ifndef GRAPHITE_IALIGNMENTMANAGER_H
#define GRAPHITE_IALIGNMENTMANAGER_H

#include "IAlignmentList.h"
#include "core/region/Region.h"

#include <boost/noncopyable.hpp>

#include <memory>

namespace graphite
{
	class Sample;
	class IAlignmentManager : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentManager > SharedPtr;
		IAlignmentManager() {}
		virtual ~IAlignmentManager() {}

        virtual IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr) = 0;
		virtual void releaseResources() = 0;
		virtual void processMappingStatistics() = 0;
		virtual std::vector< std::shared_ptr< Sample > > getSamplePtrs() = 0;
	};
}

#endif //GRAPHITE_IALIGNMENTMANAGER_H
