#ifndef GWIZ_IALIGNMENTREADERMANAGER_H
#define GWIZ_IALIGNMENTREADERMANAGER_H

#include "IAlignmentReader.h"

#include <boost/noncopyable.hpp>

#include <memory>

namespace gwiz
{
	class IAlignmentReaderManager : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentReaderManager > SharedPtr;
		IAlignmentReaderManager() {}
		virtual ~IAlignmentReaderManager() {}

		virtual IAlignmentReader::SharedPtr generateAlignmentReader() = 0;
	};
}

#endif //GWIZ_IALIGNMENTREADERMANAGER_H
