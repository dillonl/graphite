#ifndef GWIZ_ADJUDICATOR_VARIANTCONTIG_H
#define GWIZ_ADJUDICATOR_VARIANTCONTIG_H

#include <chrono>
#include <tuple>

#include <boost/noncopyable.hpp>
#include <boost/graph/depth_first_search.hpp>

#include "core/reference/IReference.h"
#include "vg/graph/VariantGraph.h"
#include "vg/graph/ReferenceNode.h"

namespace gwiz
{
namespace adjudicator
{

	class VariantContig : boost::noncopyable
	{
	public:
		typedef std::shared_ptr< VariantContig > SharedPtr;
		typedef std::tuple< std::list< vg::VariantGraph::VariantVertexDescriptor >, std::string > ContigTuple;
		typedef std::shared_ptr< ContigTuple > ContigTuplePtr;

		VariantContig(uint32_t padding, vg::VariantGraph::GraphPtr graphPtr, vg::VariantGraph::VariantVertexDescriptor startVertex, vg::VariantGraph::VariantVertexDescriptor endVertex);

		~VariantContig();

		void buildVariantContig();

		void printAllPaths()
		{
			for_each(m_contigs.begin(), m_contigs.end(), [](const ContigTuplePtr& contig){
					std::cout << std::get< 1 >(*contig) << std::endl;
				});
		}

		std::list< ContigTuplePtr > getContigs()
		{
			return m_contigs;
		}

		std::string getReferenceID();
		position getStartPosition();
		position getEndPosition();

	protected:
		void addAllPaths(std::list< vg::VariantGraph::VariantVertexDescriptor > vertexList, std::string variantSequence, vg::VariantGraph::VariantVertexDescriptor currentVertex, vg::VariantGraph::VariantVertexDescriptor endVertex);

		uint32_t m_padding;
        vg::VariantGraph::GraphPtr m_graph_ptr;
        vg::VariantGraph::VariantVertexDescriptor m_start_vertex;
        vg::VariantGraph::VariantVertexDescriptor m_end_vertex;
		std::list< ContigTuplePtr > m_contigs;
	};

} // namespace adjudicator
} // namespace gwiz

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_H
