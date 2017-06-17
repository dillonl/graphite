#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H

#include "core/alignment/IAlignmentManager.h"
#include "core/variant/IVariantManager.h"
#include "core/adjudicator/IAdjudicator.h"

#include "core/graph/GSSWGraph.h"

#include <queue>
#include <memory>
#include <mutex>

#include <unordered_map>

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
		void buildGraphs(Region::SharedPtr region, uint32_t readLength);

        std::unordered_map< std::string, std::string > getHeaderSequenceMap ();     // Return a map containing the headers and sequences for the graph paths.
        std::vector< std::string > getGraphPathHeaders ();
        std::vector< std::string > getGraphPathSequences ();
        std::vector< int > getGraphPathLengths ();

	private:
		void constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, Region::SharedPtr regionPtr, uint32_t readLength);

		std::vector< GSSWGraph::SharedPtr > m_gssw_graphs;
		std::mutex m_gssw_graph_mutex;

		IReference::SharedPtr m_reference_ptr;
		IVariantManager::SharedPtr m_variant_manager_ptr;
		IAlignmentManager::SharedPtr m_alignment_manager_ptr;
		IAdjudicator::SharedPtr m_adjudicator_ptr;

        std::unordered_map< std::string, std::string > m_header_sequence_map;  // Header sequence map for fasta.
        std::vector< std::string > m_graph_path_headers;    // Graph path headers for the resulting fasta file.
        std::vector< std::string > m_graph_path_sequences;  // Graph path sequences for the resulting fasta file.
        std::vector< int > m_graph_path_lengths;
        //std::vector< int > m_graph_path_offsets;            // Get graph path offsets to adjust the position in the resulting SAM file.
	};
}

#endif //GRAPHITE_GSSW_GSSWGRAPHMANAGER_H
