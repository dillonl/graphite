#include "core/variants/VariantList.h"
#include "GraphManager.h"

#include <queue>

namespace gwiz
{
namespace gssw
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignmentReaderManager::SharedPtr alignmentReaderManager, size_t padding) :
		m_reference_ptr(referencePtr),
		m_variant_list_ptr(variantListPtr),
		m_alignment_reader_manager(alignmentReaderManager),
		m_padding(padding)
	{
	}

	void GraphManager::buildGraphs()
	{
		auto variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
		position previousVariantEndPosition = 0;
		Variant::SharedPtr variantPtr;
		Variant::SharedPtr firstVariantInList;
		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
			if (firstVariantInList == NULL)
			{
				firstVariantInList = variantPtr;
			}
			int32_t variantPositionDifference = (variantPtr->getPosition() - previousVariantEndPosition);
			// if there is a gap between the variants that is larger than the padding then generate a graph and push it in the queue
			if ((previousVariantEndPosition > 0) && (this->m_padding < variantPositionDifference))
			{
				auto alignmentReaderPtr = this->m_alignment_reader_manager->generateAlignmentReader(); // create alignment reader

				// create region for alignmentReader
				std::string startPosition = std::to_string(firstVariantInList->getPosition() - this->m_padding);
				std::string endPosition = std::to_string(previousVariantEndPosition);
				auto region = std::make_shared< Region >(this->m_reference_ptr->getRegion()->getReferenceID() + ":" + startPosition + "-" + endPosition);

				alignmentReaderPtr->setRegion(region); // set alignmentreader's region
				auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantListPtr, alignmentReaderPtr);
				this->m_gssw_graphs.push(gsswGraph);
				variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
				firstVariantInList = variantPtr; // set the first variant
			}
			variantListPtr->addVariant(variantPtr);
			previousVariantEndPosition = variantPtr->getPosition() + variantPtr->getRef()[0].size();
		}
	}
}
}
