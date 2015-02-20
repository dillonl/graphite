#include "VariantContig.h"

#include <algorithm>
#include <deque>
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
		std::list< vg::VariantGraph::VariantVertexDescriptor > vertexList;
		std::string variantSequence;
		addAllPaths(vertexList, variantSequence, this->m_start_vertex, this->m_end_vertex);
	}

	void VariantContig::addAllPaths(std::list< vg::VariantGraph::VariantVertexDescriptor > vertexList, std::string variantSequence, vg::VariantGraph::VariantVertexDescriptor currentVertex, vg::VariantGraph::VariantVertexDescriptor endVertex)
	{
		bool finalNode = (currentVertex == endVertex);
		vertexList.push_front(currentVertex);
		auto currentNode = (*this->m_graph_ptr)[currentVertex];
		size_t length = currentNode->getLength();
		size_t offset = 0;
		if (variantSequence.empty()) // if it's the first node
		{
			offset = (this->m_region->getStartPosition() - currentNode->getPosition());
		}
		else if (finalNode) // if it's the last node
		{
			length = (this->m_region->getEndPosition() - currentNode->getPosition());
		}
		variantSequence += std::string(currentNode->getSequence() + offset, length);
		if (finalNode)
		{
			auto contig = std::make_tuple(vertexList, variantSequence);
			this->m_contigs.push_front(contig);
		}
		else
		{
			boost::graph_traits< vg::VariantGraph::Graph >::out_edge_iterator vi, vi_end;
			boost::tie(vi, vi_end) = boost::out_edges(currentVertex, *this->m_graph_ptr);
			for (; vi != vi_end; ++vi)
			{
				addAllPaths(vertexList, variantSequence, boost::target(*vi, *this->m_graph_ptr), endVertex);
			}
		}
	}

} // namespace adjudicator
} // namespace gwiz
