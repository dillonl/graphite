#ifndef GWIZ_GSSW_GSSWGRAPH_H
#define GWIZ_GSSW_GSSWGRAPH_H

#include <deque>

#include "gssw/gssw.h"

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"

namespace gwiz
{
namespace gssw
{

	class GSSWGraph : public IGraph< void >
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;

		GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;

	protected:
		std::deque< GSSWGraphPtr > m_gssw_contigs;
	};

}
}

#endif //GWIZ_GSSW_GSSWGRAPH_H
