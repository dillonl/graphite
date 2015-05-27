#ifndef GWIZ_IALIGNMENTREADER_H
#define GWIZ_IALIGNMENTREADER_H

#include <boost/noncopyable.hpp>

#include <memory>

#include "core/util/Types.h"
#include "core/region/Region.h"

#include "IAlignment.h"

namespace gwiz
{
	class IAlignmentReader : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentReader > SharedPtr;

	    IAlignmentReader()	{	}

		virtual ~IAlignmentReader()	{	}

		virtual void init() {}
		virtual void releaseReader() {}
		virtual void setRegion(Region::SharedPtr region) { this->m_region = region; }
		virtual Region::SharedPtr getRegion() { return this->m_region; }
		virtual size_t getAverageReadLength() = 0;
		virtual bool getNextAlignment(IAlignment::SharedPtr& alignment) = 0;
		virtual size_t getReadCount() = 0;

	protected:
		Region::SharedPtr m_region;

	};
}

#endif //GWIZ_IALIGNMENTREADER_H
