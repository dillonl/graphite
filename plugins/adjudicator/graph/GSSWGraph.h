#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPH_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPH_H

#include <tuple>
#include <deque>
#include <map>

#include "gssw/gssw.h"

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/allele/Allele.h"
#include "core/genotyper/IGenotyperVariant.h"
#include "core/genotyper/GenotyperAllele.hpp"

namespace graphite
{
namespace adjudicator
{

	class GSSWGraph : public IGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;
		typedef std::shared_ptr< gssw_graph_mapping > GSSWGraphMappingPtr;

		GSSWGraph(graphite::IReference::SharedPtr referencePtr, graphite::IVariantList::SharedPtr variantListPtr, position startPosition, size_t graphSize, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;
		GSSWGraphMappingPtr traceBackAlignment(IAlignment::SharedPtr alignmentPtr);
		IVariant::SharedPtr getVariantFromNodeID(const uint32_t nodeID);
		void recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr);
		gssw_graph* getGSSWGraph() { return this->m_graph_ptr; }
		int32_t getMatchValue() { return m_match; }
		IAllele::SharedPtr getAllelePtrFromNodeID(uint32_t id);
	protected:
		std::vector< gssw_node* > addAlternateVertices(const std::vector< gssw_node* >& altAndRefVertices, IVariant::SharedPtr variantPtr);
		gssw_node* addReference(position position, IAllele::SharedPtr refAllelePtr, std::vector< gssw_node* > altAndRefVertices);

		std::deque< GSSWGraphPtr > m_gssw_contigs;
		int32_t m_match;
		int32_t m_mismatch;
		int32_t m_gap_open;
		int32_t m_gap_extension;
		int8_t* m_nt_table;
		int8_t* m_mat;
		size_t m_graph_size;
		position m_start_position;
		gssw_graph* m_graph_ptr;
		uint32_t m_next_id = 0;
		std::map< uint32_t, GenotyperAllele::SharedPtr > m_genotyper_map;
		std::map< uint32_t, std::tuple< INode::SharedPtr, uint32_t, std::vector< IAlignment::SharedPtr > > > m_variant_counter;
        std::list< IGenotyperVariant > m_genotype_variants;
		std::map< uint32_t, IVariant::SharedPtr > m_variants_map;

	private:
		void graphConstructed();
		IVariantList::SharedPtr m_variant_list_ptr;
		std::vector< IAllele::SharedPtr > m_reference_fragments; // contains the reference fragments so they are deleted when the graph is deleted
		std::unordered_map< uint32_t, IAllele::SharedPtr > m_node_id_to_allele_ptrs;

		gssw_node* gssw_node_create_alt(const uint32_t position,
										const char* referenceSeq,
										const uint32_t referenceLength,
										IAllele::SharedPtr allelePtr,
										bool isReference,
										const int8_t* nt_table,
										const int8_t* score_matrix)
		{
			static std::mutex lock;
			std::lock_guard< std::mutex > guard(lock);
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			n->ref_len = referenceLength;
			n->ref_seq = (char*)referenceSeq;
			n->position = position;
			// if this node is reference then the id is even otherwise it is odd
			if (isReference)
			{
				this->m_next_id = (this->m_next_id % 2 == 0) ? this->m_next_id + 2 : this->m_next_id + 1;
				std::cout << "ref: " << std::string(allelePtr->getSequence(), allelePtr->getLength()) << std::endl;
			}
			else
			{
				std::cout << "alt: " << std::string(allelePtr->getSequence(), allelePtr->getLength()) << std::endl;
				this->m_next_id = (this->m_next_id % 2 != 0) ? this->m_next_id + 2 : this->m_next_id + 1;
			}
			n->id = this->m_next_id;
			m_node_id_to_allele_ptrs.emplace(n->id, allelePtr);
			n->len = allelePtr->getLength();
			n->seq = (char*)allelePtr->getSequence();
			n->data = (void*)allelePtr.get();
			n->num = gssw_create_num(n->seq, n->len, nt_table);
			n->count_prev = 0;
			n->count_next = 0;
			n->alignment = NULL;
			n->cigar = NULL;
			return n;
		}

	};

}
}

#endif //GRAPHITE_GSSW_GSSWGRAPH_H
