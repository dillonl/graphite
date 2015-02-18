#include "VariantContig.h"

#include <map>

namespace gwiz
{
namespace adjudicator
{

	VariantContig::VariantContig(Region::SharedPtr region, vg::VariantGraph::GraphPtr graphPtr, vg::VariantGraph::VariantVertexDescriptor startVertex, vg::VariantGraph::VariantVertexDescriptor endVertex) :
		m_region(region),
		m_graph_ptr(graphPtr),
		m_start_vertex(startVertex),
		m_end_vertex(endVertex)
	{
	}

	VariantContig::~VariantContig()
	{
	}

	void VariantContig::buildVariantContig()
	{
		this->m_contigs.clear();

		VariantContigDFSVisitor dfsVisitor;
		// boost::depth_first_search(*this->m_graph_ptr, boost::visitor(dfsVisitor).root_vertex( this->m_start_vertex), boost::get(boost::vertex_color_t(), *this->m_graph_ptr), Terminator());
		boost::depth_first_search(*this->m_graph_ptr, boost::visitor(dfsVisitor).root_vertex( this->m_start_vertex));
	}

} // namespace adjudicator
} // namespace gwiz
