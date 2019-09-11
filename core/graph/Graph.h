#ifndef GRAPHITE_GRAPH_H
#define GRAPHITE_GRAPH_H

#include "core/util/Noncopyable.hpp"
#include "core/util/GraphPrinter.h"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/vcf/Variant.h"

#include "Node.h"

#include "api/BamAlignment.h"

#include "gssw.h"

#include <memory>

namespace graphite
{
	class Graph : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Graph > SharedPtr;
		Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing, bool printGraph);
		Graph(FastaReference::SharedPtr fastaReferencePtr, Region::SharedPtr regionPtr);
		~Graph();

		std::vector< Region::SharedPtr > getRegionPtrs();
		void adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent);
        std::vector< std::vector< Node::SharedPtr > > generateAllPaths();
		Region::SharedPtr getGraphRegion();
		std::string getReferenceSequence();
		std::vector< std::string > getAllPathsAsStrings();

		void printGraphVisOutput();
		void clearResources();

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
		void processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent);
		void processTraceback2(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent);
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
		/* std::unordered_map< Node::SharedPtr, std::vector< std::string > m_paths_from_node; */
	};
}

#endif //GRAPHITE_GRAPH_H
