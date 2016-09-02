#ifndef GRAPHITE_IALIGNMENTREADER_HPP
#define GRAPHITE_IALIGNMENTREADER_HPP

#include "IAlignmentList.h"
#include "core/region/Region.h"
#include "core/util/Noncopyable.hpp"

#include <memory>

namespace graphite
{
	class IAlignmentReader : private Noncopyable
	{
	public:
		IAlignmentReader() {}
		virtual ~IAlignmentReader() {}

		virtual std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads) = 0;
	};
}

#endif //GRAPHITE_IALIGNMENTREADER_HPP
