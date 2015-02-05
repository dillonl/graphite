#ifndef GWIZ_VG_VARIANT_GRAPH_H
#define GWIZ_VG_VARIANT_GRAPH_H

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

			VariantGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr);
			virtual ~VariantGraph();

			void printGraph(const char* path);
			virtual void constructGraph() override;

			Graph::vertex_descriptor getVertexAtPosition(position referencePosition);

		protected:
			virtual bool getNextCompoundVariant(Variant::SharedPtr& variant);
			virtual Variant::SharedPtr buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< Variant::SharedPtr >& variants);

		    Graph::vertex_descriptor addReference(std::vector< Graph::vertex_descriptor >& altAndRefVertices, INode::SharedPtr referenceNode);
			std::vector< Graph::vertex_descriptor > addVariantVertices(std::vector< Graph::vertex_descriptor > altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize);
			/* position addVariant(Variant::SharedPtr variantPtr, std::vector< INode::SharedPtr >& endVertices, std::vector< INode::SharedPtr >& variantVertices); */
			/* position addReference(IReference::SharedPtr referencePtr, ); */

			GraphPtr m_graph_ptr;

			bool m_next_variant_init;
			/* typedef boost::graph_traits< Graph >::vertex_descriptor INodeType; */

			std::mutex m_graph_mutex;

			// This is used in conjunction with getCompoundNode.
			// This stores the next variant so we don't have to
			// get the variant twice
			Variant::SharedPtr m_next_variant;
			VariantParser< const char* > m_vcf_parser;

		private:



			struct OurVertexPropertyWriter {
			OurVertexPropertyWriter(Graph& g_) : g(g_) {}
				template <class Vertex>
				void operator() (std::ostream& out, Vertex u) {
					/* out << "[label=" << g[u]->nodeSeq << "]"; */
					out << "[label=\"" << g[u]->nodeSeq << " " << g[u]->seq_position << "\"]";
				}

				Graph& g;
			};
		};
	}
}

#endif // GWIZ_VG_VARIANT_GRAPH_H
