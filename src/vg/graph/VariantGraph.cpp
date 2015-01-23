#include "VariantGraph.h"

#include "ReferenceNode.h"
#include "SNPNode.h"

namespace gwiz
{
	namespace vg
	{

		VariantGraph::VariantGraph(IReference::SharedPtr referencePtr, IVariantReader::SharedPtr variantReaderPtr) :
			IGraph(referencePtr, variantReaderPtr),
			m_graph_ptr(std::make_shared< Graph >())
		{
			constructGraph();
		}

		VariantGraph::~VariantGraph()
		{
		}

		void VariantGraph::constructGraph()
		{
			// INode::SharedPtr refNode = std::make_shared< ReferenceNode >(referencePtr, )
			std::string test = "asdf";
			INode::SharedPtr node1 = std::make_shared< SNPNode >(test);
			// INode::SharedPtr node2 = std::make_shared< SNPNode >("asdf");
			// INode::SharedPtr node = std::dynamic_pointer_cast< INode >( std::make_shared< SNPNode >() );
			boost::add_vertex(node1, *m_graph_ptr);
			auto jimmer = std::dynamic_pointer_cast< SNPNode >( (*m_graph_ptr)[0] );
			std::cout << "Graph: " << jimmer->m_test << std::endl;
			// INodeType u = boost::add_vertex(node1, *m_graph_ptr);

			/*
			Graph::vertex_descriptor v0 = boost::add_vertex(*m_graph_ptr);
			Graph::vertex_descriptor v1 = boost::add_vertex(*m_graph_ptr);
			boost::add_edge(v0, v1, *m_graph_ptr);
			*/

		}

	}// end namespace vg

}// end namespace gwiz
