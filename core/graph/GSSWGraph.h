#ifndef GRAPHITE_ADJUDICATOR_GSSWGRAPH_H
#define GRAPHITE_ADJUDICATOR_GSSWGRAPH_H

#include <tuple>
#include <deque>
#include <map>
#include <mutex>

#include "gssw.h"

#include "core/graph/IGraph.h"
#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/allele/Allele.h"
#include "core/util/ThreadPool.hpp"

namespace graphite
{
	class GSSWGraphContainer
	{
	public:
	GSSWGraphContainer(int8_t* NTtable, int8_t* mat, gssw_graph* graphPtr) :
		nt_table(NTtable), mat(mat), graph_ptr(graphPtr)
		{
			lock.unlock();
		}

		~GSSWGraphContainer()
		{
			gssw_graph_destroy(this->graph_ptr);
			free(this->nt_table);
			free(this->mat);
		}

		int8_t* nt_table;
		int8_t* mat;
		gssw_graph* graph_ptr;
		std::mutex lock;
	};

	class GSSWGraph : public IGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraph > SharedPtr;
		typedef std::shared_ptr< gssw_graph > GSSWGraphPtr;
		typedef std::shared_ptr< gssw_graph_mapping > GSSWGraphMappingPtr;

		GSSWGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, Region::SharedPtr regionPtr, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue, uint32_t numGraphCopies);
		virtual ~GSSWGraph();

		virtual void constructGraph() override;
		GSSWGraphMappingPtr traceBackAlignment(IAlignment::SharedPtr alignmentPtr, std::shared_ptr< GSSWGraphContainer > graphContainer);
		/* GSSWGraphMappingPtr traceBackAlignment(IAlignment::SharedPtr alignmentPtr); */
		IVariant::SharedPtr getVariantFromNodeID(const uint32_t nodeID);
		void recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr);
		gssw_graph* getGSSWGraph() { return this->m_graph_ptr; }
		int32_t getMatchValue() { return m_match; }
		IAllele::SharedPtr getAllelePtrFromNodeID(uint32_t id);
		size_t getTotalGraphLength() { return m_total_graph_length; }
		std::string getSkipped() { return (m_skipped) ? "skipped" : "not skipped"; }

		position getStartPosition() override { this->m_region_ptr->getStartPosition(); }
		position getEndPosition() override {  this->m_region_ptr->getEndPosition(); }
		std::shared_ptr< GSSWGraphContainer > getGraphContainer();

        // Fasta getter fxns.
        std::vector< std::string > getGraphPathHeaders();
        std::vector< std::string > getGraphPathSequences();

	protected:

		void generateGraphCopies();
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
		uint32_t m_num_graph_copies;
		std::map< uint32_t, std::tuple< INode::SharedPtr, uint32_t, std::vector< IAlignment::SharedPtr > > > m_variant_counter;
		std::map< uint32_t, IVariant::SharedPtr > m_variants_map;
		std::vector< std::shared_ptr< GSSWGraphContainer > > m_graph_container_ptrs;

		size_t m_total_graph_length;
		bool m_skipped;
        
        // Fasta Headers and Sequences.
        std::vector< std::string > m_graph_path_headers;
        std::vector< std::string > m_graph_path_sequences;

		gssw_node* gssw_node_create_alt(const uint32_t position,
										const char* referenceSeq,
										const uint32_t referenceLength,
										IAllele::SharedPtr allelePtr,
										bool isReference,
										const int8_t* nt_table,
										const int8_t* score_matrix)
		{
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			n->ref_len = referenceLength;
			n->ref_seq = (char*)referenceSeq;
			n->position = position;
			// if this node is reference then the id is even otherwise it is odd
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
			n->seq = (char*)allelePtr->getSequence();
			n->data = (void*)allelePtr.get();
			n->num = gssw_create_num(n->seq, n->len, nt_table);
			n->count_prev = 0;
			n->count_next = 0;
			n->alignment = NULL;
			n->cigar = NULL;
			return n;
		}

		gssw_node* gssw_node_copy(gssw_node* node, int8_t* nt_table)
		{
			gssw_node* n = (gssw_node*)calloc(1, sizeof(gssw_node));
			n->ref_len = node->ref_len;
			n->ref_seq = node->ref_seq;
			n->position = node->position;
			n->id = node->id;
			n->len = node->len;
			n->seq = node->seq;
			n->data = node->data;
			n->num = gssw_create_num(n->seq, n->len, nt_table);
			n->count_prev = 0;
			n->count_next = 0;
			n->alignment = NULL;
			n->cigar = NULL;
			return n;
		}

	private:
		void graphConstructed();
		IVariantList::SharedPtr m_variant_list_ptr;
		std::vector< IAllele::SharedPtr > m_reference_fragments; // contains the reference fragments so they are deleted when the graph is deleted
		std::unordered_map< uint32_t, IAllele::SharedPtr > m_node_id_to_allele_ptrs;

		std::mutex m_traceback_lock;
		std::condition_variable m_condition;
		std::queue< std::shared_ptr< GSSWGraphContainer > > m_graph_container_ptrs_queue;

	};

}

#endif //GRAPHITE_GSSW_GSSWGRAPH_H
