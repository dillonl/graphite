#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPHMANAGER_H

#include "core/alignment/IAlignmentManager.h"
#include "core/variant/IVariantManager.h"
#include "core/adjudicator/IAdjudicator.h"

#include "core/graph/GSSWGraph.h"
#include "core/graph/ReferenceGraph.h"
#include "core/alignment/BamAlignment.h"

#include <queue>
#include <memory>
#include <mutex>

#include <unordered_map>

namespace graphite
{
    // May want to turn this into a struct since it only houses data.
    /**
     * This class stores the length and variantType of a node.
     */
    class NodeInfo
    {
    public:
		typedef std::shared_ptr< NodeInfo > SharedPtr;

        enum VariantType { REF, ALT };

        NodeInfo (int32_t length, VariantType variantType) :
            m_length(length),
            m_variant_type(variantType)
        {
        }

        ~NodeInfo ()
        {
        }

        int32_t getLength () { return m_length; }
        VariantType getVariantType () { return m_variant_type; }

    private:
        int32_t m_length;
        VariantType m_variant_type;
    };

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
		void buildGraphs(Region::SharedPtr region, uint32_t readLength, bool isIGVOutput);

        void registerNodeInfo(uint32_t nodeID, int length, NodeInfo::VariantType variantType);

        std::unordered_map< std::string, std::string > getHeaderSequenceMap ();     // Return a map containing the headers and sequences for the graph paths.
        std::unordered_map< uint32_t, NodeInfo::SharedPtr > getNodeInfoMap ();      // Return a map with key/value pairs of nodeIDs/NodeInfo.

	private:
		void constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, Region::SharedPtr regionPtr, uint32_t readLength, bool isIGVOutput);
        void adjudicateGraph (GSSWGraphContainer::SharedPtr gsswGraphContainer, GSSWGraphContainer::SharedPtr refGraphContainer, GSSWGraph::SharedPtr gsswGraphPtr, ReferenceGraph::SharedPtr referenceGraphPtr, IAlignment::SharedPtr alignmentPtr, position variantPosition, bool isIGVOutput);

		std::vector< GSSWGraph::SharedPtr > m_gssw_graphs;
		std::mutex m_gssw_graph_mutex;

		IReference::SharedPtr m_reference_ptr;
		IVariantManager::SharedPtr m_variant_manager_ptr;
		IAlignmentManager::SharedPtr m_alignment_manager_ptr;
		IAdjudicator::SharedPtr m_adjudicator_ptr;

        std::unordered_map< std::string, std::string > m_header_sequence_map;  // Header sequence map for fasta.
        std::unordered_map< uint32_t, NodeInfo::SharedPtr> m_node_info_map;     // Used to create bed entries.
        std::vector< graphite::BamAlignment::SharedPtr > m_bam_alignment_ptrs;    // Vector containing both the original and updated BamAlignments.
        std::vector< std::string > m_sam_entries;
	};
}

#endif //GRAPHITE_GSSW_GSSWGRAPHMANAGER_H
