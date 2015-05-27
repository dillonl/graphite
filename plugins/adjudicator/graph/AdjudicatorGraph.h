#ifndef GWIZ_ADJUDICATOR_ADJUDICATORGRAPH
#define GWIZ_ADJUDICATOR_ADJUDICATORGRAPH

#include <chrono>
#include <thread>
#include <queue>

#include "vg/graph/ReferenceNode.h"
#include "vg/graph/VariantGraph.h"

#include "core/util/ThreadPool.hpp"

#include "VariantContig.h"

#include "bamtools/src/api/BamAlignment.h"

namespace gwiz
{
namespace adjudicator
{

	class AdjudicatorGraph : public vg::VariantGraph
	{
	public:
		typedef std::shared_ptr< AdjudicatorGraph > SharedPtr;

		AdjudicatorGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr);
		virtual ~AdjudicatorGraph();

		bool MarkVariantsWithAlignment(const std::shared_ptr< BamTools::BamAlignment > alignmentPtr);

	protected:
		virtual inline vg::VariantGraph::VariantVertexDescriptor addVariantNode(INode::SharedPtr variantNodePtr) override
		{
			auto vertex = boost::add_vertex(variantNodePtr, *m_graph_ptr);
			/* this->m_reference_vertices.push_back(vertex); */
			return vertex;
		}

		inline vg::VariantGraph::VariantVertexDescriptor addReferenceNodeAtVariantPosition(vg::ReferenceNode::SharedPtr referenceNodePtr) override
		{
			auto vertex = boost::add_vertex(referenceNodePtr, *m_graph_ptr);
			this->m_reference_vertices.push_back(vertex);
			return vertex;
		}

		inline vg::VariantGraph::VariantVertexDescriptor addReferenceNode(vg::ReferenceNode::SharedPtr referenceNodePtr) override
		{
			auto vertex = boost::add_vertex(referenceNodePtr, *m_graph_ptr);
			this->m_reference_vertices.push_back(vertex);

			if (vertex == 0)
			{
				this->m_contig_start_vertex = vertex;
			}
			else if (referenceNodePtr->getLength() >= this->m_contig_padding)
			{
				auto contig = std::make_shared< VariantContig >(this->m_contig_padding, this->m_graph_ptr, this->m_contig_start_vertex, vertex);
				this->m_unprocessed_contigs.push(contig);
				this->m_contig_start_vertex = vertex;
			}
			return vertex;
		}

		inline void graphConstructed() override
		{
			processAllContigs();
		}

		void processAllContigs()
		{
			while (!this->m_unprocessed_contigs.empty())
			{
				auto contig = this->m_unprocessed_contigs.front();
				this->m_unprocessed_contigs.pop();
				contig->buildVariantContig();
				/* auto contigBuildFunct = boost::bind(&VariantContig::buildVariantContig, contig); */
				/* ThreadPool::Instance()->postJob(contigBuildFunct); */
				m_processed_contigs.push(contig);
			}
			std::cout << "finished 1" << std::endl;
			/* ThreadPool::Instance()->joinAll(); */
			std::cout << "finished 2" << std::endl;
		}

		vg::VariantGraph::VariantVertexDescriptor getReferenceVertexContainsPosition(position pos);

	public:
		std::queue< VariantContig::SharedPtr > getContigs() { return m_unprocessed_contigs; }
	private:
		const uint32_t m_contig_padding = 5; // this is the size before/after a variant that is the cut-off for variants. There must be at least m_contig_padding between contigs
		vg::VariantGraph::VariantVertexDescriptor m_contig_start_vertex;
		std::vector< vg::VariantGraph::VariantVertexDescriptor > m_reference_vertices;

		std::queue< VariantContig::SharedPtr > m_unprocessed_contigs;
		std::queue< VariantContig::SharedPtr > m_processed_contigs;

		size_t counter = 0;

	};
} // end namespace adjudicator
} // end namespace gwiz

#endif //GWIZ_ADJUDICATOR_ADJUDICATORGRAPH
