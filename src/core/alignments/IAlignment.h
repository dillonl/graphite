#ifndef GWIZ_IALIGNMENT_H
#define GWIZ_IALIGNMENT_H

#include <boost/noncopyable.hpp>

#include <memory>

#include "core/utils/Types.h"

namespace gwiz
{
	class IAlignment : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignment > SharedPtr;

	    IAlignment() {}
		virtual ~IAlignment() {}

		virtual const char* getSequence() = 0;
		virtual const position getPosition() = 0;
		virtual const size_t getLength() = 0;
		virtual const std::string getID() { return ""; }
		virtual const bool isFirstMate() { return false;}
		virtual const bool isMapped() { return false; }
		virtual const bool isReverseStrand() { return false; }
		virtual const uint16_t getOriginalMapQuality() { return 0; }

	};
}

#endif //GWIZ_IALIGNMENT_H
