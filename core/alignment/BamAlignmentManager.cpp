#include "BamAlignmentManager.h"
#include "BamAlignmentReader.h"
#include "AlignmentList.h"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <cstdlib>

namespace graphite
{
	BamAlignmentManager::BamAlignmentManager(SampleManager::SharedPtr sampleManagerPtr, Region::SharedPtr regionPtr, bool excludeDuplicateReads) :
		m_sample_manager_ptr(sampleManagerPtr),
		m_loaded(false),
        m_exclude_duplicate_reads(excludeDuplicateReads)
	{
		m_region_ptr = regionPtr;
    }

	BamAlignmentManager::BamAlignmentManager(SampleManager::SharedPtr sampleManagerPtr, Region::SharedPtr regionPtr, AlignmentReaderManager< BamAlignmentReader >::SharedPtr alignmentReaderManagerPtr, bool excludeDuplicateReads) :
		m_sample_manager_ptr(sampleManagerPtr),
		m_loaded(false),
        m_exclude_duplicate_reads(excludeDuplicateReads),
		m_alignment_reader_manager(alignmentReaderManagerPtr)
	{
		m_region_ptr = regionPtr;
    }

	BamAlignmentManager::~BamAlignmentManager()
	{
	}

	void BamAlignmentManager::loadAlignments(IVariantManager::SharedPtr variantManagerPtr)
	{
		ThreadPool::Instance()->start();
		if (this->m_loaded) { return; }
		std::unordered_set< std::string > usedPaths;
		for (auto samplePtr : this->m_sample_manager_ptr->getSamplePtrs())
		{
			if (usedPaths.find(samplePtr->getPath()) != usedPaths.end()) { continue; }
			usedPaths.emplace(samplePtr->getPath());

			std::string bamPath = samplePtr->getPath();

			std::vector< Region::SharedPtr > regionPtrs;
			auto variantListPtr = variantManagerPtr->getVariantsInRegion(this->m_region_ptr);
			IVariant::SharedPtr variantPtr;
			while (variantListPtr->getNextVariant(variantPtr))
			{
				auto variantRegionPtrs = variantPtr->getRegions();
				regionPtrs.insert(regionPtrs.end(), variantRegionPtrs.begin(), variantRegionPtrs.end());
			}

			std::deque< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;

			std::vector< BamAlignmentReader::SharedPtr > bamAlignmentReaders;
			std::unordered_set< std::string > regionStrings;
			position bamLastPosition = 0;
			bool lastPositionSet = false;
			for (auto regionPtr : regionPtrs)
			{
				auto bamAlignmentReaderPtr = m_alignment_reader_manager->getReader(bamPath);
				if (lastPositionSet)
				{
					bamLastPosition = BamAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
					lastPositionSet = true;
				}
				auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr, m_sample_manager_ptr, this->m_exclude_duplicate_reads);
				auto future = ThreadPool::Instance()->enqueue(funct);
				futureFunctions.push_back(future);
				bamAlignmentReaders.push_back(bamAlignmentReaderPtr);
			}

			{
				std::lock_guard< std::mutex > lockGuard(m_alignment_ptrs_lock);
				this->m_name_alignment_ptr_map_ptr->clear();
			}
			while (!futureFunctions.empty())
			{
				auto futureFunct = futureFunctions.front();
				futureFunctions.pop_front();
				if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
				{
					std::lock_guard< std::mutex > lockGuard(m_alignment_ptrs_lock);
					auto loadedAlignmentPtrs = futureFunct->get();
					for (auto& alignmentPtr : loadedAlignmentPtrs)
					{
						if (m_name_alignment_ptr_map_ptr->find(alignmentPtr->getID()) == m_name_alignment_ptr_map_ptr->end())
						{
							m_name_alignment_ptr_map_ptr->emplace(alignmentPtr->getID(), alignmentPtr);
							this->m_alignment_ptrs.push_back(alignmentPtr);
						}
					}
					continue;
				}
				else
				{
					futureFunctions.emplace_back(futureFunct);
				}
			}
			this->m_loaded = true;
		}

		std::sort(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(),
				  [] (const IAlignment::SharedPtr& a, const IAlignment::SharedPtr& b)
				  {
					  return a->getPosition() < b->getPosition();
				  });
	}

	void BamAlignmentManager::asyncLoadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		if (this->m_loaded) { return; }
		std::unordered_set< std::string > usedPaths;
		for (auto samplePtr : this->m_sample_manager_ptr->getSamplePtrs())
		{
			if (usedPaths.find(samplePtr->getPath()) != usedPaths.end()) { continue; }
			usedPaths.emplace(samplePtr->getPath());
			this->m_loading_thread_ptrs.emplace_back(std::make_shared< std::thread >(&BamAlignmentManager::loadBam, this, samplePtr->getPath(), variantManagerPtr, variantPadding));
		}

	}

	void BamAlignmentManager::loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		ThreadPool::Instance()->start();
		// auto regionPtrs = getRegionsContainingVariantsWithPadding(variantManagerPtr, variantPadding);
		std::vector< Region::SharedPtr > regionPtrs;
		auto variantListPtr = variantManagerPtr->getVariantsInRegion(this->m_region_ptr);
		IVariant::SharedPtr variantPtr;
		while (variantListPtr->getNextVariant(variantPtr))
		{
			auto variantRegionPtrs = variantPtr->getRegions();
			regionPtrs.insert(regionPtrs.end(), variantRegionPtrs.begin(), variantRegionPtrs.end());
		}

		std::deque< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;

		std::vector< BamAlignmentReader::SharedPtr > bamAlignmentReaders;
		std::unordered_set< std::string > regionStrings;
		for (auto regionPtr : regionPtrs)
		{
			auto bamAlignmentReaderPtr = m_alignment_reader_manager->getReader(bamPath);
			position bamLastPosition = BamAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr, m_sample_manager_ptr, this->m_exclude_duplicate_reads);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			bamAlignmentReaders.push_back(bamAlignmentReaderPtr);
		}

		{
			std::lock_guard< std::mutex > lockGuard(m_alignment_ptrs_lock);
			this->m_name_alignment_ptr_map_ptr->clear();
		}
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
			{
				std::lock_guard< std::mutex > lockGuard(m_alignment_ptrs_lock);
				auto loadedAlignmentPtrs = futureFunct->get();
				for (auto& alignmentPtr : loadedAlignmentPtrs)
				{
					if (m_name_alignment_ptr_map_ptr->find(alignmentPtr->getID()) == m_name_alignment_ptr_map_ptr->end())
					{
						m_name_alignment_ptr_map_ptr->emplace(alignmentPtr->getID(), alignmentPtr);
						this->m_alignment_ptrs.push_back(alignmentPtr);
					}
				}
				continue;
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}
		this->m_loaded = true;
	}

	SampleManager::SharedPtr BamAlignmentManager::getSamplePtrs() { return m_sample_manager_ptr; }

	void BamAlignmentManager::waitForAlignmentsToLoad()
	{
		for (auto threadPtr : this->m_loading_thread_ptrs)
		{
			threadPtr->join();
		}
		std::sort(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(),
				  [] (const IAlignment::SharedPtr& a, const IAlignment::SharedPtr& b)
				  {
					  return a->getPosition() < b->getPosition();
				  });
	}

	std::vector< Sample::SharedPtr > BamAlignmentManager::GetSamplePtrs(std::vector< std::string >& bamPaths)
	{
		std::vector< Sample::SharedPtr > samplePtrs;
		for (auto bamPath : bamPaths)
		{
			auto tmpSamplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
			samplePtrs.insert(samplePtrs.end(), tmpSamplePtrs.begin(), tmpSamplePtrs.end());
		}
		return samplePtrs;
	}

	uint32_t BamAlignmentManager::GetReadLength(std::vector< std::string >& bamPaths)
	{
		uint32_t readLength = 0;
		for (auto bamPath : bamPaths)
		{
			auto tmpReadLength = graphite::BamAlignmentReader::GetReadLength(bamPath);
			readLength = (readLength < tmpReadLength) ? tmpReadLength : readLength;
		}
		return readLength;
	}

	void BamAlignmentManager::releaseResources()
	{
	}

	/*
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
	*/
}
