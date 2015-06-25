#include "BamAlignmentManager.h"
#include "BamAlignmentReader.h"
#include "core/util/ThreadPool.hpp"

#include <functional>

namespace gwiz
{
	BamAlignmentManager::BamAlignmentManager(const std::string& bamPath, Region::SharedPtr regionPtr) :
		m_bam_path(bamPath),
		m_region_ptr(regionPtr)
	{
	}

	BamAlignmentManager::~BamAlignmentManager()
	{
	}

	IAlignmentList::SharedPtr BamAlignmentManager::getAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		// check if outside of the m_region_ptr and if so then throw an error
		return nullptr;
	}

	position BamAlignmentManager::getLastPositionInBam()
	{
		BamTools::BamReader bamReader;
		if (!bamReader.Open(this->m_bam_path))
		{
			throw "Unable to open bam file";
		}

		bamReader.LocateIndex();
		int refID = bamReader.GetReferenceID(this->m_region_ptr->getReferenceID());
		auto referenceData = bamReader.GetReferenceData();
		return referenceData[refID].RefLength;
	}

	void BamAlignmentManager::asyncLoadAlignments()
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		if (this->m_loaded) { return; }
		this->m_loading_thread_ptr = std::make_shared< std::thread >(&BamAlignmentManager::loadBam, this);
	}

	void BamAlignmentManager::loadBam()
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		std::unordered_map< std::string, bool > alignmentMap;
		std::vector< std::thread > bamLoadThreads;
		uint32_t positionIncrement = 100000;
		std::vector< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureAlignmentListPtrs;
		position startPosition = this->m_region_ptr->getStartPosition();
		position endPosition = startPosition + positionIncrement;
		position bamLastPosition = getLastPositionInBam();
		std::cout << "bamlastpos: " << bamLastPosition << std::endl;
		while (endPosition < bamLastPosition)// this->m_region_ptr->getEndPosition())
		{
            std::string regionString = this->m_region_ptr->getReferenceID() + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
			auto regionPtr = std::make_shared< Region >(regionString);
			auto bamAlignmentReaderPtr = std::make_shared< BamAlignmentReader >(this->m_bam_path);
			auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr);
			auto functFuture = ThreadPool::Instance()->enqueue(funct);
			futureAlignmentListPtrs.emplace_back(functFuture);
			startPosition = endPosition + 1;
			endPosition = ((endPosition + positionIncrement) > this->m_region_ptr->getEndPosition()) ? this->m_region_ptr->getEndPosition() : endPosition + positionIncrement;
		}
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

	void BamAlignmentManager::waitForAlignmentsToLoad()
	{
		this->m_loading_thread_ptr->join();
	}

	void BamAlignmentManager::releaseResources()
	{
	}
}
