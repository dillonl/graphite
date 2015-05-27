#ifndef GWIZ_IGRAPHADJUDICATOR_H
#define GWIZ_IGRAPHADJUDICATOR_H

#include <boost/noncopyable.hpp>
#include "core/variant/IVariantList.h"
#include "core/graph/IGraph.h"
#include "core/alignment/IAlignmentReader.h"

namespace gwiz
{
namespace gssw
{
	class IGraphAdjudicator : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IGraphAdjudicator > SharedPtr;
		IGraphAdjudicator() {}
		virtual ~IGraphAdjudicator() {}

		virtual IVariantList::SharedPtr adjudicateGraph(IGraph::SharedPtr graphPtr, IAlignmentReader::SharedPtr alignmentsPtr) = 0;
	};
}
}

#endif //GWIZ_GSSW_GRAPHADJUDICATOR_H
