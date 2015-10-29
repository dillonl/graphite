#ifndef GRAPHITE_IALIGNMENTREADER_HPP
#define GRAPHITE_IALIGNMENTREADER_HPP

#include "IAlignmentList.h"
#include "core/region/Region.h"

#include <boost/noncopyable.hpp>
#include <memory>

namespace graphite
{
	class IAlignmentReader : private boost::noncopyable
	{
	public:
		IAlignmentReader() {}
		virtual ~IAlignmentReader() {}

		virtual std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr) = 0;
	};
}

#endif //GRAPHITE_IALIGNMENTREADER_HPP
