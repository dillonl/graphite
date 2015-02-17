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
		// auto referenceVertexDescriptor = getReferenceVertexContainsPosition(alignmentPtr->Pointer);

		// loop over previously generated variant strings

		return false;
	}

	vg::VariantGraph::VariantVertexDescriptor AdjudicatorGraph::getReferenceVertexContainsPosition(position pos)
	{
		size_t startIndex = 0;
		size_t lastIndex = this->m_reference_vertices.size();
		while (startIndex <= lastIndex)
		{
			size_t midIndex = (startIndex + lastIndex) / 2;
			auto midPosition = (*this->m_graph_ptr)[this->m_reference_vertices[midIndex]]->getPosition();

			if (pos > midPosition)
			{
				startIndex = midIndex + 1;
			}
			else if (pos < midPosition)
			{
				lastIndex = midIndex - 1;
			}
			else
			{
				// return std::dynamic_pointer_cast< vg::ReferenceNode >((*this->m_graph_ptr)[this->m_reference_vertices[midIndex]]);
				return this->m_reference_vertices[midIndex];
			}
		}
		return 0;
	}

} // end namespace adjudicator
} // end namespace gwiz
