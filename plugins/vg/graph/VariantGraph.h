#ifndef GWIZ_VG_VARIANT_GRAPH_H
#define GWIZ_VG_VARIANT_GRAPH_H

#include "ReferenceNode.h"
#include "core/graph/IGraph.h"
#include "core/graph/INode.h"

#include <boost/graph/adjacency_list.hpp>

#include <mutex>

namespace gwiz
{
	namespace vg
	{
		class VariantGraph : public IGraph
		{
		public:
			typedef std::shared_ptr< VariantGraph > VariantGraphPtr;
			typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, INode::SharedPtr > Graph;
			typedef std::shared_ptr< Graph > GraphPtr;
			typedef vg::VariantGraph::Graph::vertex_descriptor VariantVertexDescriptor;

			VariantGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr);
			virtual ~VariantGraph();

			void printGraph(const char* path);
			virtual void constructGraph() override;

			VariantVertexDescriptor getVertexAtPosition(position referencePosition);

		protected:
			VariantVertexDescriptor addReference(std::vector< VariantVertexDescriptor >& altAndRefVertices, ReferenceNode::SharedPtr referenceNode);
			std::vector< VariantVertexDescriptor > addVariantVertices(std::vector< VariantVertexDescriptor > altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize);

			virtual VariantVertexDescriptor addVariantNode(INode::SharedPtr variantNodePtr)
			{
				return boost::add_vertex(variantNodePtr, *m_graph_ptr);
			}

			virtual VariantVertexDescriptor addReferenceNodeAtVariantPosition(ReferenceNode::SharedPtr referenceNodePtr)
			{
				return boost::add_vertex(referenceNodePtr, *m_graph_ptr);
			}

			virtual VariantVertexDescriptor addReferenceNode(ReferenceNode::SharedPtr referenceNodePtr)
			{
				return boost::add_vertex(referenceNodePtr, *m_graph_ptr);
			}

			/*
			 * Called when the graph is finished with construction.
			 */
			virtual inline void graphConstructed() {}

			GraphPtr m_graph_ptr;

			/* typedef boost::graph_traits< Graph >::vertex_descriptor INodeType; */

			std::mutex m_graph_mutex;

		private:



			struct OurVertexPropertyWriter {
			OurVertexPropertyWriter(Graph& g_) : g(g_) {}
				template <class Vertex>
				void operator() (std::ostream& out, Vertex u) {
					/* out << "[label=" << g[u]->nodeSeq << "]"; */
					out << "[label=\"" << std::string(g[u]->getSequence(), g[u]->getLength()) << " " << g[u]->getPosition() << "\"]";

				}

				Graph& g;
			};
		};
	}
}

#endif // GWIZ_VG_VARIANT_GRAPH_H
