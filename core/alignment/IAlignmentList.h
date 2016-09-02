#ifndef GRAPHITE_IALIGNMENTLIST_H
#define GRAPHITE_IALIGNMENTLIST_H

#include "IAlignment.h"
#include "core/region/Region.h"
#include "core/util/Noncopyable.hpp"

#include <memory>

namespace graphite
{
	class IAlignmentList : private Noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentList > SharedPtr;
		IAlignmentList() {}
		virtual ~IAlignmentList() {}

		virtual size_t getCount() = 0;
		virtual void sort() = 0;
		virtual bool getNextAlignment(IAlignment::SharedPtr& alignmentPtr) = 0;
		virtual void loadAlignmentSequences() = 0;
		virtual void unloadAlignmentSequences() = 0;
		virtual std::vector< IAlignment::SharedPtr > getAlignmentPtrs() = 0;
	protected:

	};
}

#endif //GRAPHITE_IALIGNMENTLIST_H
