#include "BamAlignmentManager.h"
#include "BamAlignmentReader.h"
#include "AlignmentList.h"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <cstdlib>

namespace graphite
{
	BamAlignmentManager::BamAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr) :
		m_region_ptr(regionPtr),
		m_sample_ptrs(samplePtrs),
		m_loaded(false)
	{
    }

	BamAlignmentManager::~BamAlignmentManager()
	{
	}

	IAlignmentList::SharedPtr BamAlignmentManager::getAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();

		auto lowerBound = std::lower_bound(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(), nullptr, [startPosition](const IAlignment::SharedPtr& alignmentPtr, const IAlignment::SharedPtr& ignore) {
				return startPosition > alignmentPtr->getPosition() + alignmentPtr->getLength();
			});
		auto upperBound = std::upper_bound(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(), nullptr, [endPosition](const IAlignment::SharedPtr& ignore, const IAlignment::SharedPtr& alignmentPtr) {
				return alignmentPtr->getPosition() > endPosition;
			});

		std::vector< IAlignment::SharedPtr > alignmentPtrs;
		alignmentPtrs.insert(alignmentPtrs.begin(), lowerBound, upperBound);
		return std::make_shared< AlignmentList >(alignmentPtrs);
	}

	void BamAlignmentManager::asyncLoadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		if (this->m_loaded) { return; }
		std::unordered_set< std::string > usedPaths;
		for (auto samplePtr : this->m_sample_ptrs)
		{
			if (usedPaths.find(samplePtr->getPath()) != usedPaths.end()) { continue; }
			usedPaths.emplace(samplePtr->getPath());
			this->m_loading_thread_ptrs.emplace_back(std::make_shared< std::thread >(&BamAlignmentManager::loadBam, this, samplePtr->getPath(), variantManagerPtr, variantPadding));
		}
	}

	std::vector< Region::SharedPtr > BamAlignmentManager::getRegionsContainingVariantsWithPadding(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		std::vector< Region::SharedPtr > regionPtrs;

		auto variantListPtr = variantManagerPtr->getVariantsInRegion(this->m_region_ptr);
		IVariant::SharedPtr variantPtr;
		IVariant::SharedPtr nextVariantPtr;
		Region::SharedPtr regionPtr = nullptr;

		position endPosition = 0;
		while (variantListPtr->getNextVariant(variantPtr))
		{
			position startPosition = ((variantPtr->getPosition() - (variantPadding * 2)) > 0) ? (variantPtr->getPosition() - (variantPadding * 2)) : 0;
			endPosition = (variantPtr->getPosition() + variantPtr->getMaxAlleleSize() + (variantPadding * 2));
			if (regionPtr == nullptr) // if this is the first time through
			{
				regionPtr = std::make_shared< Region >(variantPtr->getChrom(), startPosition, endPosition);
			}

			position regionEdgePosition = regionPtr->getEndPosition() + (variantPadding * 2);
			if (regionPtr->getReferenceID().compare(variantPtr->getChrom()) == 0 && regionEdgePosition < variantPtr->getPosition())
			{
				regionPtr->setEndPosition(endPosition); // set the end of the region variant_pos + max_allele_size + (padding * 2)
			}
			else // otherwise add the region to the list of regions and set the regionPtr = to nullptr so it will get set on the next pass
			{
				regionPtrs.push_back(regionPtr);
				regionPtr = std::make_shared< Region >(variantPtr->getChrom(), startPosition, endPosition);
			}
		}
		if (regionPtr != nullptr) // make sure to add the last region
		{
			regionPtrs.push_back(regionPtr);
		}
		return regionPtrs;
	}

	void BamAlignmentManager::loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		ThreadPool::Instance()->start();
		auto regionPtrs = getRegionsContainingVariantsWithPadding(variantManagerPtr, variantPadding);

		// std::deque< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;

		std::vector< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;
		/*
		uint32_t positionIncrement = 100000;
		for (auto regionPtr : regionPtrs)
		{
			position bamLastPosition = BamAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			// std::cout << "end position: " << regionPtr->getReferenceID() << " " << regionPtr->getStartPosition() <<  " " << regionPtr->getEndPosition() << std::endl;
			position lastRegionPosition = (bamLastPosition < regionPtr->getEndPosition()) ? bamLastPosition : regionPtr->getEndPosition();
			position startPosition = regionPtr->getStartPosition();
			position currentPosition = startPosition;
			do
			{
				auto bamAlignmentReaderPtr = std::make_shared< BamAlignmentReader >(bamPath);
				position endPosition = (lastRegionPosition < (currentPosition + positionIncrement - 1)) ? lastRegionPosition : (currentPosition + positionIncrement - 1);
				auto tmpRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), currentPosition, endPosition);
				auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, tmpRegionPtr);
				auto future = ThreadPool::Instance()->enqueue(funct);
				futureFunctions.push_back(future);
				currentPosition = endPosition + 1;
			}
			while (currentPosition < lastRegionPosition);
			// std::cout << "---" << std::endl;
		}
		*/
		// std::cout << "finished threading" << std::endl;

		std::vector< BamAlignmentReader::SharedPtr > bamAlignmentReaders;
		for (auto regionPtr : regionPtrs)
		{
			position bamLastPosition = BamAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			// std::cout << "end position: " << regionPtr->getReferenceID() << " " << regionPtr->getStartPosition() <<  " " << regionPtr->getEndPosition() << std::endl;
			// position lastRegionPosition = (bamLastPosition < regionPtr->getEndPosition()) ? bamLastPosition : regionPtr->getEndPosition();
			auto bamAlignmentReaderPtr = std::make_shared< BamAlignmentReader >(bamPath);
			auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			bamAlignmentReaders.push_back(bamAlignmentReaderPtr);
		}

		std::unordered_set< std::string > alignmentSet;
		for (auto futureFunct : futureFunctions)
		{
			futureFunct->wait();
			auto loadedAlignmentPtrs = futureFunct->get();
			for (auto& alignment : loadedAlignmentPtrs)
			{
				if (alignmentSet.find(alignment->getID()) == alignmentSet.end())
				{
					alignmentSet.insert(alignment->getID());
					this->m_alignment_ptrs.push_back(alignment);
				}
			}
		}

		/*
		std::unordered_set< std::string > alignmentSet;
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
			{
				auto loadedAlignmentPtrs = futureFunct->get();
				for (auto& alignment : loadedAlignmentPtrs)
				{
					if (alignmentSet.find(alignment->getID()) == alignmentSet.end())
					{
						alignmentSet.insert(alignment->getID());
						this->m_alignment_ptrs.push_back(alignment);
					}
				}
				// std::cout << "added: " << loadedAlignmentPtrs.size() << std::endl;
				continue;
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
			// futureFunct->wait();
		}
		*/
		// std::cout << "finished loading alignments: " << this->m_alignment_ptrs.size() << std::endl;
		this->m_loaded = true;
	}

	/*
	void BamAlignmentManager::loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		std::unordered_map< std::string, bool > alignmentMap;
		std::vector< std::thread > bamLoadThreads;
		uint32_t positionIncrement = 100000;
		std::vector< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureAlignmentListPtrs;
		position startPosition = this->m_region_ptr->getStartPosition();
		position endPosition = startPosition + positionIncrement;
		// the last position is needs to be calculated by the reader if we are loading the entire chromosome.
		// Otherwise use the specified end position.
		position bamLastPosition = (this->m_region_ptr->getEndPosition() == MAX_POSITION) ? BamAlignmentReader::GetLastPositionInBam(bamPath, this->m_region_ptr) : this->m_region_ptr->getEndPosition();

		do
		{
            std::string regionString = this->m_region_ptr->getReferenceID() + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
			auto regionPtr = std::make_shared< Region >(regionString);
			auto bamAlignmentReaderPtr = std::make_shared< BamAlignmentReader >(bamPath);
			auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr);
			auto functFuture = ThreadPool::Instance()->enqueue(funct);
			futureAlignmentListPtrs.emplace_back(functFuture);
			startPosition = endPosition + 1;
			endPosition = ((endPosition + positionIncrement) > this->m_region_ptr->getEndPosition()) ? this->m_region_ptr->getEndPosition() : endPosition + positionIncrement;
		} while (endPosition < bamLastPosition);// this->m_region_ptr->getEndPosition());
		uint32_t count = 0;
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		for (auto& alignmentFuturePtr : futureAlignmentListPtrs)
		{
			alignmentFuturePtr->wait();
			auto loadedAlignmentPtrs = alignmentFuturePtr->get();
			for (auto& alignment : loadedAlignmentPtrs)
			{
				if (alignmentMap.find(alignment->getID()) != alignmentMap.end()) { continue; } // make sure we don't add overlapping regions
				alignmentMap[alignment->getID()] = true;

				this->m_alignment_ptrs.emplace_back(alignment);
			}
		}
		this->m_loaded = true;
	}
	*/

	std::vector< Sample::SharedPtr > BamAlignmentManager::getSamplePtrs() { return m_sample_ptrs; }

	void BamAlignmentManager::waitForAlignmentsToLoad()
	{
		for (auto threadPtr : this->m_loading_thread_ptrs)
		{
			threadPtr->join();
		}
		// sort the alignments since they could be from different bam files
		std::sort(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(),
				  [] (const IAlignment::SharedPtr& a, const IAlignment::SharedPtr& b)
				  {
					  return a->getPosition() < b->getPosition();
				  });
	}

	void BamAlignmentManager::releaseResources()
	{
	}

	void BamAlignmentManager::processMappingStatistics()
	{
		for (auto alignmentPtr :  this->m_alignment_ptrs)
		{
			auto alignmentMappingMutexPtr = alignmentPtr->getMappingMutex();
			std::lock_guard< std::recursive_mutex > lock(*alignmentMappingMutexPtr); // make sure the alignmentMapping isn't set during this
			if (auto alignmentMappingWPtr = alignmentPtr->getMapping().lock()) // check if the alignment has already been aligned previously
			{
			}
		}
	}
}
