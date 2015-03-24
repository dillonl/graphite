#ifndef GWIZ_GSSW_GSSWGRAPH_H
#define GWIZ_GSSW_GSSWGRAPH_H

#include <tuple>
#include <deque>
#include <map>

#include "gssw/gssw.h"

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"
#include "core/genotyper/IGenotyperVariant.h"
#include "core/genotyper/GenotyperAllele.hpp"

#include "core/alignments/IAlignmentReader.h"

namespace gwiz
{
namespace gssw
{

	class GSSWGraph : public IGraph
	{
	public:
		static bool PrintStuff;
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;

		GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, IAlignmentReader::SharedPtr alignmentReader);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;

	protected:
		std::vector< gssw_node* > addAlternateVertices(std::vector< gssw_node* > altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize, IGenotyperVariant::SharedPtr genotyperVariantPtr);

		gssw_node* addReference(std::vector< gssw_node* > altAndRefVertices, gssw_node* referenceNode, IGenotyperVariant::SharedPtr genotyperPtr);
		gssw_node* addAlternateNode(INode::SharedPtr variantNodePtr, IGenotyperVariant::SharedPtr genotyperVariantPtr);

		std::deque< GSSWGraphPtr > m_gssw_contigs;
		IAlignmentReader::SharedPtr m_alignment_reader;
		int32_t m_match;
		int32_t m_mismatch;
		int32_t m_gap_open;
		int32_t m_gap_extension;
		int8_t* m_nt_table;
		int8_t* m_mat;
		gssw_graph* m_graph_ptr;
		uint32_t m_next_id = 0;
		position m_graph_start_position;
		std::map< uint32_t, GenotyperAllele::SharedPtr > m_genotyper_map;
		std::map< uint32_t, std::tuple< INode::SharedPtr, uint32_t, std::vector< IAlignment::SharedPtr > > > m_variant_counter;
        std::list< IGenotyperVariant > m_genotype_variants;

	private:
		void graphConstructed();
		std::shared_ptr< gssw_graph_mapping > traceBackAlignment(IAlignment::SharedPtr alignmentPtr);
		void recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr);

		gssw_node* gssw_node_create_alt(const GenotyperAllele::Type data,
										const char* seq,
										const size_t len,
										const int8_t* nt_table,
										const int8_t* score_matrix)
		{
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			//int32_t len = strlen(seq);
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
