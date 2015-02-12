#include "AdjudicatorGraph.h"

namespace gwiz
{
namespace adjudicator
{

	AdjudicatorGraph::AdjudicatorGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
		vg::VariantGraph(referencePtr, variantListPtr)
	{
	}
	AdjudicatorGraph::~AdjudicatorGraph()
	{
	}

	bool AdjudicatorGraph::MarkVariantsWithAlignment(const std::shared_ptr< BamTools::BamAlignment > alignmentPtr)
	{
		return false;
	}

	vg::ReferenceNode::SharedPtr AdjudicatorGraph::getReferenceVertexContainsPosition(position pos)
	{
		return 0;
	}

} // end namespace adjudicator
} // end namespace gwiz
