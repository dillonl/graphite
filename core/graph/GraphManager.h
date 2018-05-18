#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H

#include "core/alignment/AlignmentList.h"
#include "core/variant/VCFManager.h"
#include "core/adjudicator/GSSWAdjudicator.h"
#include "core/util/VisualizationToolKit.h"
#include "core/alignment/BamAlignmentReader.h"

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

		GraphManager(IReference::SharedPtr referencePtr, VCFManager::SharedPtr variantManagerPtr, std::vector< std::string>& bamPaths, SampleManager::SharedPtr sampleManagerPtr, bool unmappedOnly, bool includeDuplicateReads, GSSWAdjudicator::SharedPtr adjudicatorPtr);
		~GraphManager() {}

		/*
		 * graphSize: is the exact number of base pairs the graph should contain length-wise. Where there are
		 *            variants the smallest sized variant is used when adding up the graphSize.
		 *
		 * overlap: The number of base pairs that will overlap between graphs.
		 */
		void buildGraphs(Region::SharedPtr region, uint32_t readLength, VisualizationToolKit::SharedPtr vtkPtr);

	private:
		void buildGraph(std::vector< IVariant::SharedPtr > variantPtrs, position startPosition, position endPosition, Region::SharedPtr regionPtr, uint32_t readLength, VisualizationToolKit::SharedPtr vtkPtr);
		void constructAndAdjudicateGraph(VariantList::SharedPtr variantsListPtr, AlignmentList::SharedPtr alignmentListPtr, Region::SharedPtr regionPtr, uint32_t readLength, VisualizationToolKit::SharedPtr vtkPtr);

		std::vector< GSSWGraph::SharedPtr > m_gssw_graphs;
		std::mutex m_gssw_graph_mutex;

		IReference::SharedPtr m_reference_ptr;
		VCFManager::SharedPtr m_variant_manager_ptr;
		GSSWAdjudicator::SharedPtr m_adjudicator_ptr;
		std::vector< std::string> m_bam_paths;
		std::vector< BamAlignmentReader::SharedPtr > m_bam_readers;
		SampleManager::SharedPtr m_sample_manager;
		bool m_unmapped_only;
		bool m_include_duplicate_reads;
	};
}

#endif //GRAPHITE_GSSW_GSSWGRAPHMANAGER_H
