#include "core/variants/VariantList.h"
#include "GraphManager.h"

#include <queue>

namespace gwiz
{
namespace gssw
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignmentReader::SharedPtr alignmentReader, size_t padding) :
		m_reference_ptr(referencePtr),
		m_variant_list_ptr(variantListPtr),
		m_alignment_reader(alignmentReader),
		m_padding(padding)
	{
	}

	void GraphManager::buildGraphs()
	{
		auto variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
		position previousVariantEndPosition = 0;
		Variant::SharedPtr variantPtr;
		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
			int32_t variantPositionDifference = (variantPtr->getPosition() - previousVariantEndPosition);
			// if there is a gap between the variants that is larger than the padding then generate a graph and push it in the queue
			if ((previousVariantEndPosition > 0) && (this->m_padding < variantPositionDifference))
			{
				// auto gsswGraph = std::shared_ptr< GSSWGraph::SharedPtr >(this->m_reference_ptr, variantListPtr, );
				// m_gssw_graphs.push_back(gsswGraph);
				variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
			}
			variantListPtr->addVariant(variantPtr);
			previousVariantEndPosition = variantPtr->getPosition() + variantPtr->getRef()[0].size();
		}
	}
}
}
