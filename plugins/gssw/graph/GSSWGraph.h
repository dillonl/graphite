#ifndef GWIZ_GSSW_GSSWGRAPH_H
#define GWIZ_GSSW_GSSWGRAPH_H

#include <tuple>
#include <deque>
#include <map>

#include "gssw/gssw.h"

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/genotyper/IGenotyperVariant.h"
#include "core/genotyper/GenotyperAllele.hpp"

#include "core/alignment/IAlignmentReader.h"

namespace gwiz
{
namespace gssw
{

	class GSSWGraph : public IGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;
		typedef std::shared_ptr< gssw_graph_mapping > GSSWGraphMappingPtr;

		GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, position startPosition, size_t graphSize);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;
		GSSWGraphMappingPtr traceBackAlignment(IAlignment::SharedPtr alignmentPtr);
		Variant::SharedPtr getVariantFromNodeID(const uint32_t nodeID);
		void recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr);
		gssw_graph* getGSSWGraph() { return this->m_graph_ptr; }
		int32_t getMatchValue() { return m_match; }
	protected:
		std::vector< gssw_node* > addAlternateVertices(const std::vector< gssw_node* >& altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize, IGenotyperVariant::SharedPtr genotyperVariantPtr);

		gssw_node* addReference(std::vector< gssw_node* > altAndRefVertices, gssw_node* referenceNode, IGenotyperVariant::SharedPtr genotyperPtr);
		gssw_node* addAlternateNode(Variant::SharedPtr variantPtr, INode::SharedPtr variantNodePtr, IGenotyperVariant::SharedPtr genotyperVariantPtr, uint32_t variantReferenceSize);

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
		std::map< uint32_t, Variant::SharedPtr > m_variants_map;

	private:
		void graphConstructed();

		gssw_node* gssw_node_create_alt(const uint32_t position,
										const char* referenceSeq,
										const uint32_t referenceLength,
										const GenotyperAllele::Type data,
										const char* seq,
										const size_t len,
										const int8_t* nt_table,
										const int8_t* score_matrix)
		{
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			n->ref_len = referenceLength;
			n->ref_seq = (char*)malloc(n->ref_len + 1);
			strncpy(n->ref_seq, referenceSeq, n->ref_len); n->ref_seq[n->ref_len] = 0;
			n->position = position;
			n->id = m_next_id++;
			n->len = len;
			n->seq = (char*)malloc(len+1);
			strncpy(n->seq, seq, len); n->seq[len] = 0;
			n->data = (void*)data;
			n->num = gssw_create_num(seq, len, nt_table);
			n->count_prev = 0;
			n->count_next = 0;
			n->alignment = NULL;
			return n;
		}

	};

}
}

#endif //GWIZ_GSSW_GSSWGRAPH_H
