#include "GraphManager.h"
#include "ReferenceGraph.h"
#include "core/alignment/AlignmentReporter.h"
#include "core/util/ThreadPool.hpp"
#include "core/variant/VariantList.h"
#include "core/mapping/MappingManager.h"
#include "core/alignment/SampleManager.hpp"
#include "core/mapping/GSSWMapping.h"

#include "core/alignment/BamAlignmentManager.h"

#include <queue>


#include <functional>

namespace graphite
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
		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		while (currentPosition < endPosition)
		{
			auto endGraphPosition = (currentPosition + graphSize > endPosition) ? endPosition : (currentPosition + graphSize);
			auto graphRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition) + "-" + std::to_string(endGraphPosition)));
			auto variantsListPtr = this->m_variant_manager_ptr->getVariantsInRegion(graphRegion);
			if (variantsListPtr->getCount() > 0) // if we have variants, then process them
			{
				auto alignmentRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition + alignmentPadding) + "-" + std::to_string(endGraphPosition - alignmentPadding)));
				auto alignmentListPtr = this->m_alignment_manager_ptr->getAlignmentsInRegion(alignmentRegion);
				auto alignmentPtrsCount = alignmentListPtr->getCount();

				if (alignmentPtrsCount > 0)
				{
					uint32_t cutoff = ((((endGraphPosition - alignmentPadding) - (currentPosition + alignmentPadding)) / 100) * 60) * SampleManager::Instance()->getSampleCount();
					// if (alignmentPtrsCount > cutoff)
					if (false)
					{
						std::vector< IAlignment::SharedPtr > alignmentPtrs;
						IAlignment::SharedPtr alignmentPtr;
						while (alignmentListPtr->getNextAlignment(alignmentPtr))
						{
							alignmentPtrs.emplace_back(alignmentPtr);
						}
						std::random_shuffle(alignmentPtrs.begin(), alignmentPtrs.end());
						alignmentPtrs.resize(cutoff);
						auto alp = std::make_shared< AlignmentList >(alignmentPtrs);
						std::function< void() >  funct = std::bind(&GraphManager::constructAndAdjudicateGraph, this, variantsListPtr, alp, currentPosition, graphSize);
						auto future = ThreadPool::Instance()->enqueue(funct);
						futureFunctions.push_back(future);
					}
					else
					{
						std::function< void() >  funct = std::bind(&GraphManager::constructAndAdjudicateGraph, this, variantsListPtr, alignmentListPtr, currentPosition, graphSize);
						auto future = ThreadPool::Instance()->enqueue(funct);
						futureFunctions.push_back(future);
					}
				}
			}
			currentPosition += graphSize - overlap;
		}
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) != std::future_status::ready)
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}
	}

	void GraphManager::constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, position startPosition, size_t graphSize)
	{
		auto gsswGraphPtr = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, startPosition, graphSize, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue());
		gsswGraphPtr->constructGraph();

		auto referenceGraphPtr = std::make_shared< ReferenceGraph >(this->m_reference_ptr, variantsListPtr, startPosition, graphSize, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue());
		referenceGraphPtr->constructGraph();
		/*
		{
			std::lock_guard< std::mutex > lock(this->m_gssw_graph_mutex);
			this->m_gssw_graphs.emplace_back(gsswGraphPtr);
		}
		*/

		int count = 0;
		IAlignment::SharedPtr alignmentPtr;
		while (alignmentListPtr->getNextAlignment(alignmentPtr))
		{
			auto referenceMappingPtr = std::make_shared< GSSWMapping >(referenceGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
			// this->m_adjudicator_ptr->adjudicateMapping(gsswMappingPtr);

			auto referenceSWScore = referenceMappingPtr->getMappingScore();
			uint32_t referenceSWPercent = ((referenceSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);
			auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
			this->m_adjudicator_ptr->adjudicateMapping(gsswMappingPtr, referenceSWPercent);

			/*
			{
				static std::mutex m;
				std::lock_guard< std::mutex > l(m);
				gsswMappingPtr->printMapping();
			}
			*/
			MappingManager::Instance()->registerMapping(gsswMappingPtr);

		}
	}

}
