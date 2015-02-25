#include "GSSWGraph.h"


namespace gwiz
{
namespace gssw
{
	GSSWGraph::GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr) :
		IGraph(referencePtr, variantListPtr)
	{
	}

	GSSWGraph::~GSSWGraph()
	{
	}

	void GSSWGraph::constructGraph()
	{
		position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
		size_t referenceOffset = 0;
		size_t referenceSize;
		Variant::SharedPtr variantPtr;
		while (getNextCompoundVariant(variantPtr))
		{
		}
	}
}
}
