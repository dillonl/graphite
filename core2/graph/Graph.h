#ifndef GRAPHITE_GRAPH_H
#define GRAPHITE_GRAPH_H

#include "core2/util/Noncopyable.hpp"
#include "core2/region/Region.h"
#include "core2/reference/FastaReference.h"
#include "core2/vcf/Variant.h"

#include "api/BamAlignment.h"

#include "core2/util/Noncopyable.hpp"
#include "Node.h"

#include <memory>

namespace graphite
{
	class Graph : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Graph > SharedPtr;
		Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing);
		~Graph();

		std::vector< Region::SharedPtr > getRegionPtrs();
        void adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr);

	private:
		void generateGraph();
		void getGraphReference(std::string& sequence, Region::SharedPtr regionPtr);
        void generateReferenceGraphNode(Node::SharedPtr firstNodePtr, Node::SharedPtr lastNodePtr, const std::string& referenceSequence, Region::SharedPtr regionPtr);
		void addVariantsToGraph(Node::SharedPtr firstNodePtr);
		Node::SharedPtr condenseGraph(Node::SharedPtr lastNodePtr);

		FastaReference::SharedPtr m_fasta_reference_ptr;
		std::vector< Variant::SharedPtr > m_variant_ptrs;
		std::vector< Region::SharedPtr > m_graph_regions;
		uint32_t m_graph_spacing;
	};
}

#endif //GRAPHITE_GRAPH_H
