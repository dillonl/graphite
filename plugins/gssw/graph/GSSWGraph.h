#ifndef GWIZ_GSSW_GSSWGRAPH_H
#define GWIZ_GSSW_GSSWGRAPH_H

#include <deque>

#include "gssw/gssw.h"

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"

#include "core/alignments/IAlignmentReader.h"

namespace gwiz
{
namespace gssw
{

	class GSSWGraph : public IGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;

		GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, IAlignmentReader::SharedPtr alignmentReader);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;
		void constructAndAlign();

	protected:
		virtual inline void graphConstructed() {}
		/* virtual void  */
		std::deque< GSSWGraphPtr > m_gssw_contigs;

		IAlignmentReader::SharedPtr m_alignment_reader;
		int32_t m_match;
		int32_t m_mismatch;
		int32_t m_gap_open;
		int32_t m_gap_extension;
		int8_t* m_nt_table;
		int8_t* m_mat;
		gssw_graph* m_graph_ptr;

	private:
		bool recenterGraph(IAlignment::SharedPtr alignment);
		void constructGraph(std::list< Variant::SharedPtr > variants, size_t referencePadding);

		gssw_node* gssw_node_create(const char* seq,
									const size_t len,
									const int8_t* nt_table,
									const int8_t* score_matrix) {
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			/* int32_t len = strlen(seq); */
			/* n->id = id; */
			n->len = len;
			n->seq = (char*)malloc(len+1);
			strncpy(n->seq, seq, len); n->seq[len] = 0;
			/* n->data = data; */
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
