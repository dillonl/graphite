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
		vg::ReferenceNode::SharedPtr getReferenceVertexContainsPosition(position pos);

	};
} // end namespace adjudicator
} // end namespace gwiz

#endif //GWIZ_ADJUDICATOR_ADJUDICATORGRAPH
