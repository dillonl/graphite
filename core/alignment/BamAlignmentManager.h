#ifndef GRAPHITE_BAMALIGNMENTMANAGER_H
#define GRAPHITE_BAMALIGNMENTMANAGER_H

#include "IAlignmentManager.h"
#include "core/region/Region.h"
#include "core/variant/IVariantList.h"
#include "core/variant/IVariantManager.h"
#include "Sample.hpp"

#include <mutex>
#include <thread>
#include <atomic>
#include <memory>
#include <future>

namespace graphite
{
	class BamAlignmentManager : public IAlignmentManager
	{
	public:
		typedef std::shared_ptr< BamAlignmentManager > SharedPtr;
		BamAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr, bool excludeDuplicateReads = false);
		~BamAlignmentManager();

		void asyncLoadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);
		void waitForAlignmentsToLoad();
		void releaseResources() override;
		void processMappingStatistics() override;
		std::vector< Sample::SharedPtr > getSamplePtrs() override;
	private:
		void loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);

		std::mutex m_loaded_mutex;
		bool m_loaded;
        bool m_exclude_duplicate_reads;
		std::string m_bam_path;
		std::vector< std::shared_ptr< std::thread > > m_loading_thread_ptrs;
		std::mutex m_alignment_ptrs_lock;
		std::vector< Sample::SharedPtr > m_sample_ptrs;
	};
}

#endif //GRAPHITE_BAMALIGNMENTMANAGER_H
