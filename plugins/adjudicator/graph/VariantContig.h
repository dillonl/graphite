#ifndef GWIZ_ADJUDICATOR_VARIANTCONTIG_H
#define GWIZ_ADJUDICATOR_VARIANTCONTIG_H

#include <boost/noncopyable.hpp>

#include <boost/graph/depth_first_search.hpp>

#include "core/reference/IReference.h"
#include "vg/graph/VariantGraph.h"
#include "vg/graph/ReferenceNode.h"

namespace gwiz
{
namespace adjudicator
{

	class VariantContigDFSVisitor : public boost::default_dfs_visitor
	{
    public:
		template < typename Vertex, typename Graph >
		void discover_vertex(const Vertex& v, const Graph& g) const
		{

		}
	};

	class Terminator
	{
    public:
		template < typename Vertex, typename Graph >
		bool operator()(const Vertex& v, const Graph& g) const
		{
			return false;
		}
	};

	class VariantContig : boost::noncopyable
	{
	public:
VariantContig(Region::SharedPtr region, vg::VariantGraph::GraphPtr graphPtr, vg::VariantGraph::VariantVertexDescriptor startVertex, vg::VariantGraph::VariantVertexDescriptor endVertex);
		~VariantContig();

		void buildVariantContig();

	private:
		Region::SharedPtr m_region;
        vg::VariantGraph::GraphPtr m_graph_ptr;
        vg::VariantGraph::VariantVertexDescriptor m_start_vertex;
        vg::VariantGraph::VariantVertexDescriptor m_end_vertex;
		std::vector< std::string > m_contigs;
	};

} // namespace adjudicator
} // namespace gwiz

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_H
