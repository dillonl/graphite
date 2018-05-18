#ifndef GRAPHITE_GRAPH_REFERENCEGRAPH_H
#define GRAPHITE_GRAPH_REFERENCEGRAPH_H

#include "GSSWGraph.h"
#include "core/reference/IReference.h"
#include "core/variant/VariantList.h"
#include "core/allele/Allele.h"

#include "gssw.h"


namespace graphite
{
	class ReferenceGraph : public GSSWGraph
	{
	public:
		typedef std::shared_ptr< ReferenceGraph > SharedPtr;

	    ReferenceGraph(IReference::SharedPtr referencePtr, VariantList::SharedPtr variantListPtr, Region::SharedPtr regionPtr, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue);
		virtual ~ReferenceGraph();

		void constructGraph() override;

	};
}

#endif //GRAPHITE_GR_GSSWGRAPHMANAGER_H
