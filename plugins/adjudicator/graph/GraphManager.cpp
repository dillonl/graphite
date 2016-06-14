#include "GraphManager.h"
#include "AlignmentReporter.h"
#include "core/util/ThreadPool.hpp"
#include "core/genotyper/IGenotyper.h"
#include "core/variant/VariantList.h"
#include "core/mapping/MappingManager.h"
#include "GSSWMapping.h"

#include "core/alignment/BamAlignmentManager.h"

#include <queue>


#include <boost/bind.hpp>

namespace graphite
{
namespace adjudicator
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantManager::SharedPtr variantManagerPtr, IAlignmentManager::SharedPtr alignmentManagerPtr, IAdjudicator::SharedPtr adjudicatorPtr) :
		m_reference_ptr(referencePtr),
		m_variant_manager_ptr(variantManagerPtr),
		m_alignment_manager_ptr(alignmentManagerPtr),
		m_adjudicator_ptr(adjudicatorPtr)
	{
	}

	void GraphManager::buildGraphs(Region::SharedPtr regionPtr, size_t graphSize, size_t overlap, size_t alignmentPadding)
	{
		std::string referenceID = regionPtr->getReferenceID();
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();
		position currentPosition = startPosition;

		while (currentPosition < endPosition)
		{
			auto endGraphPosition = (currentPosition + graphSize > endPosition) ? endPosition : (currentPosition + graphSize);
			auto graphRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition) + "-" + std::to_string(endGraphPosition)));
			auto variantsListPtr = this->m_variant_manager_ptr->getVariantsInRegion(graphRegion);
			if (variantsListPtr->getCount() > 0) // if we have variants, then process them
			{
				auto alignmentRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition + alignmentPadding) + "-" + std::to_string(endGraphPosition - alignmentPadding)));
				auto alignmentListPtr = this->m_alignment_manager_ptr->getAlignmentsInRegion(alignmentRegion);
				if (alignmentListPtr->getCount() > 0)
				{
					auto funct = std::bind(&GraphManager::constructAndAdjudicateGraph, this, variantsListPtr, alignmentListPtr, currentPosition, graphSize);
					ThreadPool::Instance()->enqueue(funct);
				}
			}
			currentPosition += graphSize - overlap;
		}
		ThreadPool::Instance()->joinAll();
	}

	void GraphManager::constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, position startPosition, size_t graphSize)
	{
		auto gsswGraphPtr = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, startPosition, graphSize, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue());
		gsswGraphPtr->constructGraph();
		{
			std::lock_guard< std::mutex > lock(this->m_gssw_graph_mutex);
			this->m_gssw_graphs.emplace_back(gsswGraphPtr);
		}

		// std::cout << "graph manager -- load sequences: " << startPosition << "-" << (startPosition+graphSize) << std::endl;
		// alignmentListPtr->loadAlignmentSequences();
		IAlignment::SharedPtr alignmentPtr;
		while (alignmentListPtr->getNextAlignment(alignmentPtr))
		{
			auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
			MappingManager::Instance()->registerMapping(gsswMappingPtr);
		}
		// alignmentListPtr->unloadAlignmentSequences();
	}

}
}
