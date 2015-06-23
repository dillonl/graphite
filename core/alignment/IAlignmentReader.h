#ifndef GWIZ_IALIGNMENTREADER_HPP
#define GWIZ_IALIGNMENTREADER_HPP

#include "IAlignmentList.h"
#include "core/region/Region.h"

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class IAlignmentReader : private boost::noncopyable
	{
	public:
		IAlignmentReader() {}
		virtual ~IAlignmentReader() {}

		virtual IAlignmentList::SharedPtr loadAlignmentsInRegion(Region::SharedPtr regionPtr) = 0;
	};
}

#endif //GWIZ_IALIGNMENTREADER_HPP
