#ifndef GWIZ_IALIGNMENTREADER_H
#define GWIZ_IALIGNMENTREADER_H

#include <boost/noncopyable.hpp>

#include <memory>

#include "core/utils/Types.h"
#include "core/region/Region.h"

#include "IAlignment.h"

namespace gwiz
{
	class IAlignmentReader : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentReader > SharedPtr;

	    IAlignmentReader() {}
		virtual ~IAlignmentReader() {}

		virtual void setRegion(Region::SharedPtr region) = 0;
		virtual size_t getAverageReadLength() = 0;
		virtual bool getNextAlignment(IAlignment::SharedPtr& alignment) = 0;

	};
}

#endif //GWIZ_IALIGNMENTREADER_H
