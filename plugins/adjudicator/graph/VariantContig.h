#ifndef GWIZ_ADJUDICATOR_VARIANTCONTIG_H
#define GWIZ_ADJUDICATOR_VARIANTCONTIG_H

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
VariantContig(Region::SharedPtr region, vg::VariantGraph::GraphPtr graphPtr, vg::VariantGraph::VariantVertexDescriptor startVertex, vg::VariantGraph::VariantVertexDescriptor endVertex);
		~VariantContig();

		void buildVariantContig();

		void printAllPaths()
		{
			for_each(m_contigs.begin(), m_contigs.end(), [](const std::tuple< std::list< vg::VariantGraph::VariantVertexDescriptor >, std::string >& contig){
					std::cout << std::get< 1 >(contig) << std::endl;
				});
		}

	protected:
		void addAllPaths(std::list< vg::VariantGraph::VariantVertexDescriptor > vertexList, std::string variantSequence, vg::VariantGraph::VariantVertexDescriptor currentVertex, vg::VariantGraph::VariantVertexDescriptor endVertex);

		Region::SharedPtr m_region;
        vg::VariantGraph::GraphPtr m_graph_ptr;
        vg::VariantGraph::VariantVertexDescriptor m_start_vertex;
        vg::VariantGraph::VariantVertexDescriptor m_end_vertex;
		std::list< std::tuple< std::list< vg::VariantGraph::VariantVertexDescriptor >, std::string > > m_contigs;
	};

} // namespace adjudicator
} // namespace gwiz

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_H
