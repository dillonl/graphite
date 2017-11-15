#include "GraphManager.h"
#include "ReferenceGraph.h"
#include "core/alignment/AlignmentReporter.h"
#include "core/util/ThreadPool.hpp"
#include "core/variant/VariantList.h"
#include "core/mapping/MappingManager.h"
#include "core/mapping/GSSWMapping.h"

#include "core/alignment/BamAlignmentManager.h"

#include <queue>
#include <algorithm>
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

	void GraphManager::buildGraphs(Region::SharedPtr regionPtr, uint32_t readLength)
	{
		auto variantsListPtr = this->m_variant_manager_ptr->getVariantsInRegion(regionPtr);
		if (variantsListPtr->getCount() == 0) // if we don't have variants or alignments in the region, then return
		{
			return;
		}

		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		// loop through variants and build and adjudicate graphs
		IVariant::SharedPtr variantPtr = nullptr;
		while (variantsListPtr->getNextVariant(variantPtr))
		{
			if (variantPtr->shouldSkip()) { continue; }
			position startPosition = 0;
			position endPosition = 0;
			auto variantRegionPtrs = variantPtr->getRegions();
			if (variantRegionPtrs.size() > 0)
			{
				startPosition = variantPtr->getRegions()[0]->getStartPosition();
				endPosition = variantPtr->getRegions()[variantPtr->getRegions().size() - 1]->getEndPosition();
			}
			std::vector< IVariant::SharedPtr > variantPtrs;
			IVariant::SharedPtr nextVariantPtr = nullptr;
			variantPtrs.emplace_back(variantPtr);
			while (variantsListPtr->peekNextVariant(nextVariantPtr) && variantPtr->doesOverlap(nextVariantPtr, readLength) && (!variantPtr->isStructuralVariant() && !nextVariantPtr->isStructuralVariant()))
			{
				variantsListPtr->getNextVariant(nextVariantPtr);
				variantPtrs.emplace_back(nextVariantPtr);

				for (auto regionPtr : nextVariantPtr->getRegions())
				{
					auto tmpStartPosition = regionPtr->getStartPosition();
					auto tmpEndPosition = regionPtr->getEndPosition();
					startPosition = (tmpStartPosition < startPosition) ? tmpStartPosition : startPosition;
					endPosition = (tmpEndPosition > endPosition) ? tmpEndPosition : endPosition;
				}
			}
			// getting all the alignments in the variant's region
			std::vector< IAlignment::SharedPtr > alignmentPtrs;
			std::unordered_set< IAlignment* > alignmentPtrSet;
			for (auto regionVariantPtr : variantPtrs)
			{
				auto regionPtrs = regionVariantPtr->getRegions();
				for (auto regionPtr : regionPtrs)
				{
					auto regionAlignmentListPtr = this->m_alignment_manager_ptr->getAlignmentsInRegion(regionPtr);
					auto regionAlignmentPtrsCount = regionAlignmentListPtr->getCount();
					if (regionAlignmentPtrsCount > 0)
					{
						auto regionAlignmentPtrs = regionAlignmentListPtr->getAlignmentPtrs();
						// alignmentPtrs.insert(alignmentPtrs.end(), regionAlignmentPtrs.begin(), regionAlignmentPtrs.end()); // combine lists of alignments
						for (auto alignmentPtr : regionAlignmentPtrs)
						{
							if (alignmentPtrSet.find(alignmentPtr.get()) == alignmentPtrSet.end())
							{
								alignmentPtrSet.emplace(alignmentPtr.get());
								alignmentPtrs.emplace_back(alignmentPtr);
								if (alignmentPtr->getPosition() < startPosition) { startPosition = alignmentPtr->getPosition(); }
								if ((alignmentPtr->getPosition() + alignmentPtr->getLength()) > endPosition) { endPosition = alignmentPtr->getPosition() + alignmentPtr->getLength(); }
							}
						}
					}
				}
			}
			alignmentPtrSet.clear();
			// std::cout << regionPtr->getRegionString() << " " << alignmentPtrs.size() << std::endl;
			if (alignmentPtrs.size() > 0)
			{
				// find the start and end position for the graph
				auto graphAlignmentRegion = std::make_shared< Region >(regionPtr->getReferenceID(), startPosition, endPosition, Region::BASED::ONE);
				for (auto i = 1; i < alignmentPtrs.size(); ++i)
				{
					if (alignmentPtrs[i]->getPosition() < startPosition) { startPosition = alignmentPtrs[i]->getPosition(); }
					if (alignmentPtrs[i]->getPosition() + alignmentPtrs[0]->getLength() > endPosition) { endPosition = alignmentPtrs[i]->getPosition() + alignmentPtrs[0]->getLength(); }
				}
				auto variantListPtr = std::make_shared< VariantList >(variantPtrs, this->m_reference_ptr);
				auto alignmentListPtr = std::make_shared< AlignmentList >(alignmentPtrs);
				constructAndAdjudicateGraph(variantListPtr, alignmentListPtr, graphAlignmentRegion, readLength);
			}
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

	void GraphManager::constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, Region::SharedPtr regionPtr, uint32_t readLength)
	{
		// uint32_t numGraphCopies = (alignmentListPtr->getCount() < ThreadPool::Instance()->getThreadCount()) ? alignmentListPtr->getCount() : ThreadPool::Instance()->getThreadCount();  // get the min of threadcount and alignment count, this is the num of simultanious threads processing this graph
		uint32_t numGraphCopies = ThreadPool::Instance()->getThreadCount();  // get the min of threadcount and alignment count, this is the num of simultanious threads processing this graph
		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		auto gsswGraphPtr = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, regionPtr, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue(), numGraphCopies);
		gsswGraphPtr->constructGraph();

		auto referenceGraphPtr = std::make_shared< ReferenceGraph >(this->m_reference_ptr, variantsListPtr, regionPtr, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue(), numGraphCopies);
		referenceGraphPtr->constructGraph();

		/*
		static int count = 0;

		std::cout << std::string(40, '-') << std::endl;
		std::cout << regionPtr->getRegionString() << std::endl;

		for (auto variantPtr : variantsListPtr->getAllVariantPtrs())
		{
			std::cout << variantPtr->getPosition() << " " << variantPtr->getRefAllelePtr()->getSequenceString();
			for (auto altAllelePtr : variantPtr->getAltAllelePtrs())
			{
				std::cout << " " << altAllelePtr->getSequenceString();
			}
			std::cout << std::endl;
		}
		std::cout << "===PATHS===" << std::endl;

		auto paths = std::static_pointer_cast< GSSWGraph >(gsswGraphPtr)->generateAllPaths();
		// get all graph paths from Andrew
		for (auto i = 0; i < paths.size(); ++i)
		{
			std::cout << std::get< 1 >(paths[i]) << std::endl;
			std::cout << std::get< 0 >(paths[i]) << std::endl;
			std::cout << std::endl;
		}
		std::cout << "===ALIGNMENTS===" << std::endl;
		*/

		IAlignment::SharedPtr alignmentPtr;
		auto alignmentPtrs = alignmentListPtr->getAlignmentPtrs();
		while (alignmentListPtr->getNextAlignment(alignmentPtr))
		{
			/*
			static std::mutex mut;
			{
				mut.lock();
				std::cout << std::string(alignmentPtr->getSequence(), alignmentPtr->getLength()) << " " << alignmentPtr->getID().substr(0, alignmentPtr->getID().size() - 1) << std::endl;
				mut.unlock();
			}
			*/
			auto gsswGraphContainer = gsswGraphPtr->getGraphContainer();
			auto refGraphContainer = referenceGraphPtr->getGraphContainer();
			auto funct = [gsswGraphContainer, refGraphContainer, gsswGraphPtr, referenceGraphPtr, alignmentPtr, this]()
			{
				auto refTraceback = referenceGraphPtr->traceBackAlignment(alignmentPtr, refGraphContainer);
				auto referenceMappingPtr = std::make_shared< GSSWMapping >(refTraceback, alignmentPtr);
				auto referenceSWScore = referenceMappingPtr->getMappingScore();
				uint32_t referenceSWPercent = ((referenceSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);

				auto tracebackPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr, gsswGraphContainer);
				auto gsswMappingPtr = std::make_shared< GSSWMapping >(tracebackPtr, alignmentPtr);

				auto gsswSWScore = referenceMappingPtr->getMappingScore();
				uint32_t gsswSWPercent = ((gsswSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);

				if (this->m_adjudicator_ptr->adjudicateMapping(gsswMappingPtr, referenceSWPercent))
				{
					MappingManager::Instance()->registerMapping(gsswMappingPtr);
				}
		    };

			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
		}

		// wait for all the functions to complete (for all graphs to be created and adjudicated)
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

}
