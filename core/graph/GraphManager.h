#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H

#include "core/alignment/IAlignmentManager.h"
#include "core/variant/IVariantManager.h"
#include "core/adjudicator/IAdjudicator.h"

#include "core/graph/GSSWGraph.h"

#include <queue>
#include <memory>
#include <mutex>

namespace graphite
{
	class GraphManager : private Noncopyable
	{
	public:
		typedef std::shared_ptr< GraphManager > SharedPtr;

		GraphManager(IReference::SharedPtr referencePtr, IVariantManager::SharedPtr variantManagerPtr, IAlignmentManager::SharedPtr alignmentManagerPtr, IAdjudicator::SharedPtr adjudicatorPtr);
		~GraphManager() {}

		/*
		 * graphSize: is the exact number of base pairs the graph should contain length-wise. Where there are
		 *            variants the smallest sized variant is used when adding up the graphSize.
		 *
		 * overlap: The number of base pairs that will overlap between graphs.
		 */
		void buildGraphs(Region::SharedPtr region, size_t graphSize, size_t overlap, size_t alignmentPadding);

	private:
		void constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, position startPosition, size_t graphSize);
		void adjudicateList(std::vector< IAlignment::SharedPtr > alignmentPtrs, GSSWGraph::SharedPtr graphPtr);

		std::vector< GSSWGraph::SharedPtr > m_gssw_graphs;
		std::mutex m_gssw_graph_mutex;

		IReference::SharedPtr m_reference_ptr;
		IVariantManager::SharedPtr m_variant_manager_ptr;
		IAlignmentManager::SharedPtr m_alignment_manager_ptr;
		IAdjudicator::SharedPtr m_adjudicator_ptr;
	};
}

#endif //GRAPHITE_GSSW_GSSWGRAPHMANAGER_H
