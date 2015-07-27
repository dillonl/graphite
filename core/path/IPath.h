#ifndef GWIZ_IPATH_H
#define GWIZ_IPATH_H

#include <boost/noncopyable.h>

namespace gwiz
{
	class IPath : boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IPath > SharedPtr;
		IPath() {}
		~IPath() {}

		virtual std::vector< IAllele::SharedPtr > getAllelePath() = 0;
	};
}

#endif //GWIZ_IPATH_H
