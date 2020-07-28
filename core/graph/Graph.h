#ifndef GRAPHITE_GRAPH_H
#define GRAPHITE_GRAPH_H

#include "core/util/Noncopyable.hpp"
#include "core/util/GraphPrinter.h"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/vcf/Variant.h"
#include "core/alignment/Alignment.h"

#include "Node.h"

#include "gssw.h"

#include <memory>

namespace graphite
{
	class Graph
	{
	public:
		typedef std::shared_ptr< Graph > SharedPtr;
		Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing, bool printGraph);
		Graph(FastaReference::SharedPtr fastaReferencePtr, Region::SharedPtr regionPtr);
		Graph(const Graph& graph) = delete;
		Graph& operator=(const Graph& graph) = delete;
		Graph() {}
		~Graph();

		std::vector< Region::SharedPtr > getRegionPtrs();
        std::vector< std::vector< Node::SharedPtr > > generateAllPaths();
		Region::SharedPtr getGraphRegion();
		std::string getReferenceSequence();
		std::vector< std::string > getAllPathsAsStrings();

		void printGraphVisOutput();
		void clearResources();
		Graph::SharedPtr createCopy();
		std::unordered_map< uint32_t, Node::SharedPtr > getNodePtrsMap() { return m_node_ptrs_map; }
		Node::SharedPtr getFirstNode() { return m_first_node; }

		void removeNodePtr(Node* nodePtr);

	private:
		std::vector< Variant::SharedPtr > reconcileVariantSemantics(std::vector< Variant::SharedPtr >& variantPtrs);
		void addAdditionalNodeEdges();
		void compressLargeNodes();
		void setRegionPtrs();
		void generateGraph();
		void generateReferenceGraph(Region::SharedPtr regionPtr);
		void getGraphReference(std::string& sequence, Region::SharedPtr& regionPtr, std::vector< Variant::SharedPtr >& variantPtrs);
        void generateReferenceGraphNode(Node::SharedPtr& firstNodePtr, Node::SharedPtr& lastNodePtr, const std::string& referenceSequence, Region::SharedPtr regionPtr);
		void addVariantsToGraph(Node::SharedPtr firstNodePtr);
		Node::SharedPtr condenseGraph(Node::SharedPtr lastNodePtr);
		std::vector< std::string > generateAllPathsFromNodesOfLength(Node::SharedPtr nodePtr);
		void setPrefixAndSuffix(Node::SharedPtr firstNodePtr);

		bool isNodePrefixAmbiguous(std::string& nodeSequence, Node* comparatorNode, std::unordered_set< Node::SharedPtr >& nodePtrs);
		bool isNodeSuffixAmbiguous(std::string& nodeSequence, Node* comparatorNode, std::unordered_set< Node::SharedPtr >& nodePtrs);

		FastaReference::SharedPtr m_fasta_reference_ptr;
		std::vector< Variant::SharedPtr > m_variant_ptrs;
		std::vector< Region::SharedPtr > m_graph_regions;
		std::unordered_map< uint32_t, Node::SharedPtr > m_node_ptrs_map;
		uint32_t m_graph_spacing;
		Node::SharedPtr m_first_node;
		uint32_t m_score_threshold;
		std::unordered_set< Node::SharedPtr > m_all_created_nodes;
		std::unordered_set< std::string > m_aligned_read_names;
		std::mutex m_aligned_read_names_mutex;
        GraphPrinter::SharedPtr m_graph_printer_ptr;
		std::mutex m_graph_mutex;
		/* std::unordered_map< Node::SharedPtr, std::vector< std::string > m_paths_from_node; */
	};

}

#endif //GRAPHITE_GRAPH_H
