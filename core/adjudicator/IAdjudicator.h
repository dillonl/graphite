#ifndef GRAPHITE_IADJUDICATOR_H
#define GRAPHITE_IADJUDICATOR_H

#include "core/mapping/IMapping.h"
#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class IAdjudicator : private Noncopyable, public std::enable_shared_from_this< IAdjudicator >
	{
	public:
		typedef std::shared_ptr< IAdjudicator > SharedPtr;
		IAdjudicator() {}
		virtual ~IAdjudicator() {}

		virtual void adjudicateMapping(IMapping::SharedPtr mappingPtr, uint32_t referenceSWPercent) = 0;
		virtual int getMatchValue() = 0;
		virtual int getMisMatchValue() = 0;
		virtual int getGapOpenValue() = 0;
		virtual int getGapExtensionValue() = 0;
	private:
	};
}

#endif //GRAPHITE_ADJUDICATOR_H
