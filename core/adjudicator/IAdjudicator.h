#ifndef GWIZ_IADJUDICATOR_H
#define GWIZ_IADJUDICATOR_H

#include "core/mapping/IMapping.h"

namespace gwiz
{
	class IAdjudicator : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAdjudicator > SharedPtr;
		IAdjudicator() {}
		virtual ~IAdjudicator() {}

		virtual void adjudicateMapping(IMapping::SharedPtr mappingPtr) = 0;
		virtual int getMatchValue() = 0;
		virtual int getMisMatchValue() = 0;
		virtual int getGapOpenValue() = 0;
		virtual int getGapExtensionValue() = 0;
	private:
	};
}

#endif //GWIZ_ADJUDICATOR_H
