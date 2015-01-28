#ifndef GWIZ_VG_VARIANT_GRAPH_H
#define GWIZ_VG_VARIANT_GRAPH_H

//#include "core/variants/IVariantIterator.h"
#include "core/graph/IGraph.h"
#include "core/graph/INode.h"

/* #include <boost/graph/graph_traits.hpp> */
#include <boost/graph/adjacency_list.hpp>

namespace gwiz
{
	namespace vg
	{
		class VariantGraph : public IGraph
		{
		public:
			typedef std::shared_ptr< VariantGraph > VariantGraphPtr;

			VariantGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr);
			~VariantGraph();

			void printGraph(std::string& path);

		protected:
			void constructGraph() override;

		private:
			typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, INode::SharedPtr > Graph;
			typedef std::shared_ptr< Graph > GraphPtr;
			/* typedef boost::graph_traits< Graph >::vertex_descriptor INodeType; */

			GraphPtr m_graph_ptr;


			struct OurVertexPropertyWriter {
			OurVertexPropertyWriter(Graph& g_) : g(g_) {}
				template <class Vertex>
				void operator() (std::ostream& out, Vertex u) {
					out << "[label=" << g[u]->nodeSeq << "]";
				}

				Graph& g;
			};
		};
	}
}

#endif // GWIZ_VG_VARIANT_GRAPH_H
