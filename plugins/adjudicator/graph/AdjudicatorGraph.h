#ifndef GWIZ_ADJUDICATOR_ADJUDICATORGRAPH
#define GWIZ_ADJUDICATOR_ADJUDICATORGRAPH

#include "vg/graph/ReferenceNode.h"
#include "vg/graph/VariantGraph.h"

#include "bamtools/src/api/BamAlignment.h"

namespace gwiz
{
namespace adjudicator
{

	class AdjudicatorGraph : public vg::VariantGraph
	{
	public:
		typedef std::shared_ptr< AdjudicatorGraph > SharedPtr;

		AdjudicatorGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr);
		virtual ~AdjudicatorGraph();

		bool MarkVariantsWithAlignment(const std::shared_ptr< BamTools::BamAlignment > alignmentPtr);

	protected:
		inline vg::VariantGraph::VariantVertexDescriptor addReferenceNode(vg::ReferenceNode::SharedPtr referenceNodePtr) override
		{

			auto vertex = boost::add_vertex(referenceNodePtr, *m_graph_ptr);
			this->m_reference_vertices.push_back(vertex);
			return vertex;
		}

		vg::VariantGraph::VariantVertexDescriptor getReferenceVertexContainsPosition(position pos);

	protected:

	private:
		std::vector< vg::VariantGraph::VariantVertexDescriptor > m_reference_vertices;

	};
} // end namespace adjudicator
} // end namespace gwiz

#endif //GWIZ_ADJUDICATOR_ADJUDICATORGRAPH
