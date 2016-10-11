#include "HTSLibAlignmentManager.h"
#include "HTSLibAlignmentReader.h"
// #include "BamAlignmentReader.h"
#include "AlignmentList.h"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <cstdlib>

namespace graphite
{
	HTSLibAlignmentManager::HTSLibAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr, AlignmentReaderManager< HTSLibAlignmentReader >::SharedPtr alignmentReaderManagerPtr, bool excludeDuplicateReads) :
		m_sample_ptrs(samplePtrs),
        m_exclude_duplicate_reads(excludeDuplicateReads),
		m_loaded(false),
		m_alignment_reader_manager(alignmentReaderManagerPtr)
	{
		m_region_ptr = regionPtr;
	}
	HTSLibAlignmentManager::~HTSLibAlignmentManager()
	{
	}

    IAlignmentList::SharedPtr HTSLibAlignmentManager::getAlignmentsInRegion(Region::SharedPtr regionPtr)
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
		return std::make_shared< AlignmentList >(alignmentPtrs,m_name_alignment_ptr_map_ptr, regionPtr);
	}

	void HTSLibAlignmentManager::releaseResources()
	{
	}

	void HTSLibAlignmentManager::processMappingStatistics()
	{
	}

	void HTSLibAlignmentManager::loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		ThreadPool::Instance()->start();
		auto regionPtrs = getRegionsContainingVariantsWithPadding(variantManagerPtr, variantPadding);

		std::deque< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;

		std::vector< HTSLibAlignmentReader::SharedPtr > alignmentReaders;
		std::unordered_set< std::string > regionStrings;

		for (auto regionPtr : regionPtrs)
		{
			// position bamLastPosition = BamAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			auto alignmentReaderPtr = m_alignment_reader_manager->getReader(bamPath);
			position bamLastPosition = alignmentReaderPtr->getRegionLastPosition(regionPtr);
			// auto alignmentReaderPtr = std::make_shared< HTSLibAlignmentReader >(bamPath);
			// alignmentReaderPtr->open();
			// auto funct = std::bind(&HTSLibAlignmentReader::loadAlignmentsInRegion, alignmentReaderPtr, regionPtr, m_alignment_reader_manager, this->m_exclude_duplicate_reads);
			auto funct = std::bind(&HTSLibAlignmentReader::loadAlignmentsInRegion, alignmentReaderPtr, regionPtr, this->m_exclude_duplicate_reads);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			alignmentReaders.push_back(alignmentReaderPtr);
		}

		// std::unordered_set< std::string > alignmentSet;
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
				// std::cout << "loadedAlignments: " << loadedAlignmentPtrs.size()<< std::endl;
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

	void HTSLibAlignmentManager::loadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		ThreadPool::Instance()->start();
		if (this->m_loaded) { return; }

		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;
		std::unordered_set< std::string > usedPaths;
		for (auto samplePtr : this->m_sample_ptrs)
		{
			if (usedPaths.find(samplePtr->getPath()) != usedPaths.end()) { continue; }
			usedPaths.emplace(samplePtr->getPath());
			auto funct = std::bind(&HTSLibAlignmentManager::loadBam, this, samplePtr->getPath(), variantManagerPtr, variantPadding);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			// this->m_loading_thread_ptrs.emplace_back(std::make_shared< std::thread >(&HTSLibAlignmentManager::loadBam, this, samplePtr->getPath(), variantManagerPtr, variantPadding));
		}

		for (auto futureFunct : futureFunctions)
		{
			futureFunct->wait();
		}
		/*
		for (auto threadPtr : this->m_loading_thread_ptrs)
		{
			threadPtr->join();
		}
		*/
		// sort the alignments since they could be from different bam files
		std::sort(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(),
				  [] (const IAlignment::SharedPtr& a, const IAlignment::SharedPtr& b)
				  {
					  return a->getPosition() < b->getPosition();
				  });
	}
}
