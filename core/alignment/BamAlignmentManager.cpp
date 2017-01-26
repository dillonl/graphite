#include "BamAlignmentManager.h"
#include "BamAlignmentReader.h"
#include "SamtoolsAlignmentReader.h"
#include "AlignmentList.h"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <cstdlib>

namespace graphite
{
	BamAlignmentManager::BamAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr, bool excludeDuplicateReads) :
		m_sample_ptrs(samplePtrs),
		m_loaded(false),
        m_exclude_duplicate_reads(excludeDuplicateReads)
	{
		m_region_ptr = regionPtr;
    }

	BamAlignmentManager::BamAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr, AlignmentReaderManager< BamAlignmentReader >::SharedPtr alignmentReaderManagerPtr, bool excludeDuplicateReads) :
		m_sample_ptrs(samplePtrs),
		m_loaded(false),
        m_exclude_duplicate_reads(excludeDuplicateReads),
		m_alignment_reader_manager(alignmentReaderManagerPtr)
	{
		m_region_ptr = regionPtr;
    }

	BamAlignmentManager::~BamAlignmentManager()
	{
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

	void BamAlignmentManager::loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		ThreadPool::Instance()->start();
		auto regionPtrs = getRegionsContainingVariantsWithPadding(variantManagerPtr, variantPadding);

		std::deque< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;

		// std::vector< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;


		// std::vector< SamtoolsAlignmentReader::SharedPtr > bamAlignmentReaders;
		std::vector< BamAlignmentReader::SharedPtr > bamAlignmentReaders;
		std::unordered_set< std::string > regionStrings;
		for (auto regionPtr : regionPtrs)
		{
			// position bamLastPosition = SamtoolsAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			auto bamAlignmentReaderPtr = m_alignment_reader_manager->getReader(bamPath);
			// std::cout << "Using bam reader: " << bamAlignmentReaderPtr->getReaderID() << std::endl;
			position bamLastPosition = BamAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			// auto bamAlignmentReaderPtr = std::make_shared< SamtoolsAlignmentReader >(bamPath);
			// auto bamAlignmentReaderPtr = std::make_shared< BamAlignmentReader >(bamPath);
			// bamAlignmentReaderPtr->open();
			// auto funct = std::bind(&SamtoolsAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr, this->m_exclude_duplicate_reads);
			auto funct = std::bind(&BamAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr, this->m_exclude_duplicate_reads);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			bamAlignmentReaders.push_back(bamAlignmentReaderPtr);
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
				// std::cout << "alignments size: " << loadedAlignmentPtrs.size() << std::endl;
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



		/*

		  for testing only
		for (auto alignmentPtr : this->m_alignment_ptrs)
		{
			std::cout << alignmentPtr->getPosition() << " " << alignmentPtr->getSequence() << std::endl;
		}
		*/
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
