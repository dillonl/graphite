#include "BamAlignmentManager.h"
#include "BamAlignmentReader.h"
#include "AlignmentList.h"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <boost/filesystem.hpp>
#include <unordered_set>

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

	void BamAlignmentManager::asyncLoadAlignments()
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		if (this->m_loaded) { return; }
		std::unordered_set< std::string > usedPaths;
		for (auto samplePtr : this->m_sample_ptrs)
		{
			if (usedPaths.find(samplePtr->getPath()) != usedPaths.end()) { continue; }
			this->m_loading_thread_ptr = std::make_shared< std::thread >(&BamAlignmentManager::loadBam, this, samplePtr->getPath());
		}
	}

	void BamAlignmentManager::loadBam(const std::string& bamPath)
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
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

	std::vector< Sample::SharedPtr > BamAlignmentManager::getSamplePtrs() { return m_sample_ptrs; }

	void BamAlignmentManager::waitForAlignmentsToLoad()
	{
		this->m_loading_thread_ptr->join();
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
