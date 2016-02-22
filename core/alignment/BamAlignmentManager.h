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
		BamAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr);
		~BamAlignmentManager();

		IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr) override;
		void asyncLoadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);
		void loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);
		void waitForAlignmentsToLoad();
		void releaseResources() override;
		void processMappingStatistics() override;
		std::vector< Sample::SharedPtr > getSamplePtrs() override;
	private:
		std::vector< Region::SharedPtr > getRegionsContainingVariantsWithPadding(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);

		std::mutex m_loaded_mutex;
		bool m_loaded;
		std::string m_bam_path;
		Region::SharedPtr m_region_ptr;
		std::vector< std::shared_ptr< std::thread > > m_loading_thread_ptrs;
		std::mutex m_alignment_ptrs_lock;
        std::vector< IAlignment::SharedPtr > m_alignment_ptrs;
		std::vector< Sample::SharedPtr > m_sample_ptrs;
	};
}

#endif //GRAPHITE_BAMALIGNMENTMANAGER_H
