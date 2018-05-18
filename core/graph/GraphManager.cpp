#include "GraphManager.h"
#include "ReferenceGraph.h"
#include "core/alignment/AlignmentReporter.h"
#include "core/variant/VariantList.h"
#include "core/mapping/GSSWMapping.h"

#include <queue>
#include <algorithm>
#include <functional>

namespace graphite
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, VCFManager::SharedPtr variantManagerPtr, std::vector< std::string>& bamPaths, SampleManager::SharedPtr sampleManagerPtr, bool unmappedOnly, bool includeDuplicateReads, GSSWAdjudicator::SharedPtr adjudicatorPtr) :
		m_reference_ptr(referencePtr),
		m_variant_manager_ptr(variantManagerPtr),
		m_bam_paths(bamPaths),
		m_adjudicator_ptr(adjudicatorPtr),
		m_sample_manager(sampleManagerPtr),
		m_unmapped_only(unmappedOnly),
		m_include_duplicate_reads(includeDuplicateReads)
	{
		for (auto bamPath : m_bam_paths)
		{
			auto bamReaderPtr = std::make_shared< BamAlignmentReader >(bamPath);
			bamReaderPtr->open();
			m_bam_readers.emplace_back(bamReaderPtr);
		}
	}

	void GraphManager::buildGraphs(Region::SharedPtr regionPtr, uint32_t readLength, VisualizationToolKit::SharedPtr vtkPtr)
	{
		auto variantsListPtr = this->m_variant_manager_ptr->getVariantsInRegion(regionPtr);
		if (variantsListPtr->getCount() == 0) // if we don't have variants or alignments in the region, then return
		{
			return;
		}

		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		// loop through variants and build and adjudicate graphs
		IVariant::SharedPtr variantPtr = nullptr;

		std::deque< std::shared_future< void > > futures;
		int numThreads = std::thread::hardware_concurrency() * 2;

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
			bool exitLoop = false;
			while (!exitLoop)
			{
				if (futures.size() < numThreads)
				{
					auto f = std::async(std::launch::async, std::bind(&GraphManager::buildGraph, this, variantPtrs, startPosition, endPosition, regionPtr, readLength, vtkPtr));
					futures.emplace_back(f.share());
					break;
				}
				else
				{
					auto f = futures.front();
					futures.pop_front();
					if (f.wait_for(std::chrono::milliseconds(1)) == std::future_status::ready)
					{
						auto fut = std::async(std::launch::async, std::bind(&GraphManager::buildGraph, this, variantPtrs, startPosition, endPosition, regionPtr, readLength, vtkPtr));
						auto futShare = fut.share();
						futures.emplace_back(futShare);
						break;
					}
					futures.emplace_back(f);
				}
			}
		}
		while (!futures.empty())
		{
			auto f = futures.front();
			futures.pop_front();
			f.wait();
		}
	}

	void GraphManager::buildGraph(std::vector< IVariant::SharedPtr > variantPtrs, position startPosition, position endPosition, Region::SharedPtr regionPtr, uint32_t readLength, VisualizationToolKit::SharedPtr vtkPtr)
	{
		// getting all the alignments in the variant's region
		std::vector< IAlignment::SharedPtr > alignmentPtrs;
		std::unordered_set< std::string > alignmentPtrSet;
		for (auto regionVariantPtr : variantPtrs)
		{
			auto regionPtrs = regionVariantPtr->getRegions();
			for (auto regionPtr : regionPtrs)
			{
				position beforeStartPosition = regionPtr->getStartPosition() - 1000;
				position afterEndPosition = regionPtr->getEndPosition() + 1000;
				auto beforeRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), beforeStartPosition, regionPtr->getStartPosition(), regionPtr->getBased());
				auto afterRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), regionPtr->getEndPosition(), afterEndPosition, regionPtr->getBased());
				auto graphAlignmentRegion1 = std::make_shared< Region >(regionPtr->getReferenceID(), startPosition, endPosition, Region::BASED::ONE);


				{
					for (BamAlignmentReader::SharedPtr bamReaderPtr : m_bam_readers)
					{
						BamAlignmentReader bamReader(bamReaderPtr->getPath());
						bamReader.open();
						auto tmpAlignmentPtrsBefore = bamReader.loadAlignmentsInRegion(beforeRegionPtr, m_sample_manager, true,  m_include_duplicate_reads);
						auto tmpAlignmentPtrs = bamReader.loadAlignmentsInRegion(regionPtr, m_sample_manager, false,  m_include_duplicate_reads);
						auto tmpAlignmentPtrsAfter = bamReader.loadAlignmentsInRegion(afterRegionPtr, m_sample_manager, true,  m_include_duplicate_reads);

						auto addAlignmentFunct = [&alignmentPtrSet, &alignmentPtrs](std::vector< IAlignment::SharedPtr >& alignmentListPtrs)
							{
								for (IAlignment::SharedPtr alignmentPtr : alignmentListPtrs)
								{
									if (alignmentPtrSet.find(alignmentPtr->getID()) == alignmentPtrSet.end())
									{
										alignmentPtrSet.emplace(alignmentPtr->getID());
										alignmentPtrs.emplace_back(alignmentPtr);
									}
								}
							};
						addAlignmentFunct(tmpAlignmentPtrsBefore);
						addAlignmentFunct(tmpAlignmentPtrs);
						addAlignmentFunct(tmpAlignmentPtrsAfter);
						bamReader.close();
					}
				}
			}

		}
		alignmentPtrSet.clear();
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
			constructAndAdjudicateGraph(variantListPtr, alignmentListPtr, graphAlignmentRegion, readLength, vtkPtr);
		}
	}

	void GraphManager::constructAndAdjudicateGraph(VariantList::SharedPtr variantsListPtr, AlignmentList::SharedPtr alignmentListPtr, Region::SharedPtr regionPtr, uint32_t readLength, VisualizationToolKit::SharedPtr vtkPtr)
	{
		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		IAlignment::SharedPtr alignmentPtr;
		auto alignmentPtrs = alignmentListPtr->getAlignmentPtrs();
		while (alignmentListPtr->getNextAlignment(alignmentPtr))
		{

			auto gsswGraphPtr = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, regionPtr, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue());
			gsswGraphPtr->constructGraph();

			auto referenceGraphPtr = std::make_shared< ReferenceGraph >(this->m_reference_ptr, variantsListPtr, regionPtr, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue());
			referenceGraphPtr->constructGraph();

			auto refTraceback = referenceGraphPtr->traceBackAlignment(alignmentPtr);

			auto referenceMappingPtr = std::make_shared< GSSWMapping >(refTraceback, alignmentPtr, referenceGraphPtr->getStartPosition());
			auto referenceSWScore = referenceMappingPtr->getMappingScore();
			uint32_t referenceSWPercent = ((referenceSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);

			auto tracebackPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr);
			auto gsswMappingPtr = std::make_shared< GSSWMapping >(tracebackPtr, alignmentPtr, gsswGraphPtr->getStartPosition());

			auto gsswSWScore = referenceMappingPtr->getMappingScore();
			uint32_t gsswSWPercent = ((gsswSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);

			if (this->m_adjudicator_ptr->adjudicateMapping(gsswMappingPtr, referenceSWPercent) && vtkPtr != nullptr)
			{
				vtkPtr->setAlignmentAndMapping(alignmentPtr, gsswGraphPtr, referenceMappingPtr, gsswMappingPtr);
			}
		}
	}

}
