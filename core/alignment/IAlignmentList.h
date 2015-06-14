#ifndef GWIZ_IALIGNMENTLIST_H
#define GWIZ_IALIGNMENTLIST_H

#include "IAlignment.h"
#include "core/region/Region.h"

#include <boost/noncopyable.hpp>

#include <memory>

namespace gwiz
{
	class IAlignmentList : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentList > SharedPtr;
		IAlignmentList() {}
		virtual ~IAlignmentList() {}

		virtual size_t getCount() = 0;
		virtual void sort() = 0;
		virtual bool getNextAlignment(IAlignment::SharedPtr& alignmentPtr) = 0;
	protected:

	};
}

#endif //GWIZ_IALIGNMENTLIST_H
