#include "GraphManager.h"
#include "AlignmentReporter.h"
#include "core/util/ThreadPool.hpp"
#include "core/genotyper/IGenotyper.h"
#include "core/variant/VariantList.h"

#include "core/alignment/BamAlignmentManager.h"

#include <queue>


#include <boost/bind.hpp>

namespace gwiz
{
namespace gssw
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantManager::SharedPtr variantManagerPtr, IAlignmentManager::SharedPtr alignmentManagerPtr, IGraphAdjudicator::SharedPtr graphAdjudicatorPtr) :
		m_reference_ptr(referencePtr),
		m_variant_manager_ptr(variantManagerPtr),
		m_alignment_manager_ptr(alignmentManagerPtr),
		m_graph_adjudicator_ptr(graphAdjudicatorPtr)
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
				// auto alignmentRegion = std::make_shared< Region >("1:12307541-12309541");
				// create region for alignmentReader
				// alignmentReaderPtr->init();
				// alignmentReaderPtr->setRegion(alignmentRegion); // set alignmentreader's region
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
		auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, startPosition, graphSize);
		gsswGraph->constructGraph();
		this->m_gssw_graphs.emplace_back(gsswGraph);
		this->m_graph_adjudicator_ptr->adjudicateGraph(gsswGraph, alignmentListPtr);
	}

}
}
