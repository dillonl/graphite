#ifndef GWIZ_GSSW_GSSWGRAPH_H
#define GWIZ_GSSW_GSSWGRAPH_H

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"

namespace gwiz
{
namespace gssw
{

	class GSSWGraph : public IGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;

		GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;

	protected:

	};

}
}

#endif //GWIZ_GSSW_GSSWGRAPH_H
