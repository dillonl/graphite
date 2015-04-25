#ifndef GWIZ_GSSW_GSSWGRAPHMANAGER_H
#define GWIZ_GSSW_GSSWGRAPHMANAGER_H

#include "core/alignments/IAlignmentReaderManager.h"
#include "core/variants/VariantList.h"
#include "core/adjudicator/IGraphAdjudicator.h"

#include "GSSWGraph.h"

#include <queue>
#include <memory>
#include <mutex>

#include <boost/noncopyable.hpp>

namespace gwiz
{
namespace gssw
{
	class GraphManager : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< GraphManager > SharedPtr;

		GraphManager(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, IAlignmentReaderManager::SharedPtr alignmentReaderManager, IGraphAdjudicator::SharedPtr graphAdjudicatorPtr);
		~GraphManager() {}

		/*
		 * graphSize: is the exact number of base pairs the graph should contain length-wise. Where there are
		 *            variants the smallest sized variant is used when adding up the graphSize.
		 *
		 * overlap: The number of base pairs that will overlap between graphs.
		 */
		IVariantList::SharedPtr buildGraphs(Region::SharedPtr region, size_t graphSize, size_t overlap, size_t alignmentPadding);

	private:
		/* void constructAndAdjudicateGraph(IVariantList::SharedPtr reportedVariants, IVariantList::SharedPtr variantsListPtr, IAlignmentReader::SharedPtr alignmentReaderPtr, position startPosition, size_t graphSize,std::mutex& reportedVariantsMutex); */
		void constructAndAdjudicateGraph(IVariantList::SharedPtr reportedVariants, IVariantList::SharedPtr variantsListPtr, IAlignmentReader::SharedPtr alignmentReaderPtr, position startPosition, size_t graphSize);

		std::queue< GSSWGraph::SharedPtr > m_gssw_graphs;
		std::mutex m_gssw_graph_mutex;

		IReference::SharedPtr m_reference_ptr;
		IVariantList::SharedPtr m_variant_list_ptr;
		IAlignmentReaderManager::SharedPtr m_alignment_reader_manager;
		IGraphAdjudicator::SharedPtr m_graph_adjudicator_ptr;
	};
}
}

#endif //GWIZ_GSSW_GSSWGRAPHMANAGER_H
