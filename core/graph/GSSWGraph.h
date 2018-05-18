#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPH_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPH_H

#include <tuple>
#include <deque>
#include <map>
#include <mutex>

#include "gssw.h"

#include "core/reference/IReference.h"
#include "core/variant/VariantList.h"
#include "core/allele/Allele.h"

namespace graphite
{

	class GSSWGraph : private Noncopyable
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;
		typedef std::shared_ptr< gssw_graph_mapping > GSSWGraphMappingPtr;

		GSSWGraph(IReference::SharedPtr referencePtr, VariantList::SharedPtr variantListPtr, Region::SharedPtr regionPtr, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue);
		virtual ~GSSWGraph();

		virtual void constructGraph();
		GSSWGraphMappingPtr traceBackAlignment(IAlignment::SharedPtr alignmentPtr);
		IVariant::SharedPtr getVariantFromNodeID(const uint32_t nodeID);
		void recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr);
		gssw_graph* getGSSWGraph() { return this->m_graph_ptr; }
		int32_t getMatchValue() { return m_match; }
		IAllele::SharedPtr getAllelePtrFromNodeID(uint32_t id);
		size_t getTotalGraphLength() { return m_total_graph_length; }
		std::string getSkipped() { return (m_skipped) ? "skipped" : "not skipped"; }

		position getStartPosition() { this->m_region_ptr->getStartPosition(); }
		position getEndPosition() {  this->m_region_ptr->getEndPosition(); }

		std::vector< std::tuple< std::string, std::string > > generateAllPaths();

	protected:

		std::vector< gssw_node* > addAlternateVertices(const std::vector< gssw_node* >& altAndRefVertices, IVariant::SharedPtr variantPtr);
		gssw_node* addReferenceVertex(position position, IAllele::SharedPtr refAllelePtr, std::vector< gssw_node* > altAndRefVertices);

		std::deque< GSSWGraphPtr > m_gssw_contigs;
		int32_t m_match;
		int32_t m_mismatch;
		int32_t m_gap_open;
		int32_t m_gap_extension;
		int8_t* m_nt_table;
		int8_t* m_mat;
		gssw_graph* m_graph_ptr;
		Region::SharedPtr m_region_ptr;
		static uint32_t s_next_id;
		static std::mutex s_lock;
		std::map< uint32_t, IVariant::SharedPtr > m_variants_map;

		size_t m_total_graph_length;
		bool m_skipped;

		gssw_node* gssw_node_create_alt(const uint32_t position,
										const char* referenceSeq,
										const uint32_t referenceLength,
										IAllele::SharedPtr allelePtr,
										bool isReference,
										const int8_t* nt_table,
										const int8_t* score_matrix)
		{
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			/* if this node is reference then the id is even otherwise it is odd */
			char* tmpSeq = (char*)malloc(allelePtr->getLength() + 1 * sizeof(char));
			memcpy(tmpSeq, allelePtr->getSequence(), allelePtr->getLength() + 1);
			{
				std::lock_guard< std::mutex > l(s_lock);
				if (isReference)
				{
					s_next_id = (s_next_id % 2 == 0) ? s_next_id + 2 : s_next_id + 1;
				}
				else
				{
					s_next_id = (s_next_id % 2 != 0) ? s_next_id + 2 : s_next_id + 1;
				}
				n->id = s_next_id;
			}
			allelePtr->setID(n->id);

			if (m_node_id_to_allele_ptrs.find(n->id) == m_node_id_to_allele_ptrs.end())
			{
				m_node_id_to_allele_ptrs.emplace(n->id, allelePtr);
			}
			n->len = allelePtr->getLength();
			n->seq = tmpSeq;
			n->data = (void*)allelePtr.get();
			allelePtr->setPosition(position);
			n->num = gssw_create_num(n->seq, n->len, nt_table);
			n->count_prev = 0;
			n->count_next = 0;
			n->alignment = NULL;
			return n;
		}

		gssw_node* gssw_node_copy(gssw_node* node, int8_t* nt_table)
		{
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			n->id = node->id;
			n->seq = node->seq;
			n->len = node->len;
			n->data = node->data;
			n->num = gssw_create_num(n->seq, n->len, nt_table);
			n->count_prev = 0;
			n->count_next = 0;
			n->alignment = NULL;
			return n;
		}

		IReference::SharedPtr m_reference_ptr;

	private:
		void graphConstructed();
		VariantList::SharedPtr m_variant_list_ptr;
		std::vector< IAllele::SharedPtr > m_reference_fragments; // contains the reference fragments so they are deleted when the graph is deleted
		std::unordered_map< uint32_t, IAllele::SharedPtr > m_node_id_to_allele_ptrs;

	};

}

#endif //GRAPHITE_GSSW_GSSWGRAPH_H
