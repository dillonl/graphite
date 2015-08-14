#include <boost/graph/graphviz.hpp>
#include <boost/algorithm/string/join.hpp>

#include "VariantGraph.h"
#include "SNPNode.h"
#include "PathTraceVisitor.h"

#include <deque>
#include <set>

namespace graphite
{
	namespace vg
	{

		VariantGraph::VariantGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
			IGraph(referencePtr, variantListPtr),
			m_graph_ptr(std::make_shared< Graph >())
		{
		}

		VariantGraph::~VariantGraph()
		{
		}

		void VariantGraph::getAllPaths(std::vector< std::string >& paths, std::vector< std::vector< INode::SharedPtr > >& nodes)
		{
			Graph graph = *this->m_graph_ptr;
			std::deque< VariantVertexDescriptor > vertices;
			vertices.emplace_front(boost::vertex(0,graph));
			std::vector< INode::SharedPtr > nodeTracker;
			std::unordered_map< VariantVertexDescriptor, std::vector< INode::SharedPtr > > vertexMap;
			while (!vertices.empty())
			{
				auto currentVertexDescriptor = vertices.front();
				vertices.pop_front();
				auto nodePtr = graph[currentVertexDescriptor];

				if (vertexMap.find(currentVertexDescriptor) != vertexMap.end())
				{
					nodeTracker = vertexMap[currentVertexDescriptor];
				}
				nodeTracker.emplace_back(nodePtr);
				Graph::out_edge_iterator e, end;
				boost::tie(e, end) = boost::out_edges(currentVertexDescriptor, graph);
				bool pathEnd = (e == end);
				for (; e != end; ++e)
				{
					auto vertexDescriptor = target(*e, graph);
					vertices.emplace_front(vertexDescriptor);
					vertexMap[vertexDescriptor] = nodeTracker;
				}
				if (pathEnd)
				{
					nodes.emplace_back(nodeTracker);
				}
			}
		}

		VariantGraph::VariantVertexDescriptor VariantGraph::getVertexAtPosition(position referencePosition)
		{
			throw "UNIMPLEMENTED";
			// size_t startIndex = 0;
			// size_t lastIndex = this->m_reference_vertices.size() - 1;
			// while (startIndex <= lastIndex)
			// {
			// 	size_t midIndex = (startIndex + lastIndex) / 2;
			// 	auto midPosition = (*this->m_graph_ptr)[this->m_reference_vertices[midIndex]]->getPosition();
			// 	if (pos > midPosition)
			// 	{
			// 		startIndex = midIndex + 1;
			// 	}
			// 	else if (pos < midPosition)
			// 	{
			// 		lastIndex = midIndex - 1;
			// 	}
			// 	else
			// 	{
			// 		return this->m_reference_vertices[midIndex];
			// 	}
			// }
			// return this->m_reference_vertices[lastIndex];
		}

		void VariantGraph::constructGraph()
		{
			position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
			size_t referenceOffset = 0;
			size_t referenceSize;
			VariantVertexDescriptor referenceVertex;
			std::vector< VariantVertexDescriptor > altAndRefVertices;
			IVariant::SharedPtr variantPtr;
			while (this->m_variant_list_ptr->getNextVariant(variantPtr))
			{
				referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
				if (referenceSize > 0) // there is reference to add
				{
					ReferenceNode::SharedPtr referenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, referenceSize);
					referenceVertex = addReference(altAndRefVertices, referenceNode);
					altAndRefVertices.clear();
					altAndRefVertices.push_back(referenceVertex);
				}
				size_t variantReferenceSize;
				altAndRefVertices = addVariantVertices(altAndRefVertices, variantPtr, variantReferenceSize);
				referenceOffset += referenceSize + variantReferenceSize;
			}
			referenceSize = this->m_reference_ptr->getRegion()->getEndPosition() - (startPosition + referenceOffset);
			if (referenceSize > 0)
			{
				ReferenceNode::SharedPtr referenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, referenceSize);
				referenceVertex = addReference(altAndRefVertices, referenceNode);
			}
			graphConstructed();
		}

		VariantGraph::VariantVertexDescriptor VariantGraph::addReference(std::vector< VariantGraph::VariantVertexDescriptor >& altAndRefVertices, ReferenceNode::SharedPtr referenceNode)
		{
			VariantVertexDescriptor referenceVertex = addReferenceNode(referenceNode);
			// add previous variant and reference Vertices to the referenceVertex
			for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
			{
				boost::add_edge((*iter), referenceVertex, *this->m_graph_ptr);
			}
			return referenceVertex;
		}

		std::vector< VariantGraph::VariantVertexDescriptor > VariantGraph::addVariantVertices(std::vector< VariantGraph::VariantVertexDescriptor > altAndRefVertices, IVariant::SharedPtr variantPtr, size_t& variantReferenceSize)
		{
			std::vector< VariantVertexDescriptor > vertices;
			for (uint32_t i = 0; i < variantPtr->getAltAllelePtrs().size(); ++i)
			{
				INode::SharedPtr variantNode = IVariantNode::BuildVariantNodes(variantPtr, i);
				vertices.push_back(addVariantNode(variantNode));
			}
			size_t referenceOffset = variantPtr->getPosition() - this->m_reference_ptr->getRegion()->getStartPosition();
			size_t refSize = variantPtr->getRefAllelePtr()->getSequenceString().size();
			ReferenceNode::SharedPtr variantReferenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, refSize);
			vertices.push_back(addReferenceNodeAtVariantPosition(variantReferenceNode));
			variantReferenceSize = refSize;

			for (auto parentVertexIter = altAndRefVertices.begin(); parentVertexIter != altAndRefVertices.end(); ++parentVertexIter)
			{
				for (auto childVertexIter = vertices.begin(); childVertexIter != vertices.end(); ++childVertexIter)
				{
					boost::add_edge((*parentVertexIter), (*childVertexIter), *this->m_graph_ptr);
				}
			}

			return vertices;
		}

		void VariantGraph::printGraph(const char* path)
		{
			std::ofstream ofs(path);
			boost::write_graphviz(ofs, *this->m_graph_ptr,OurVertexPropertyWriter(*this->m_graph_ptr));
		}

	}// end namespace vg

}// end namespace graphite
