#ifndef GWIZ_VCFMANAGER_H
#define GWIZ_VCFMANAGER_H

#include "core/region/Region.h"
#include "VCFList.h"
#include "VariantList.h"
#include "IVariantManager.h"
#include "IVariant.h"

#include <atomic>
#include <memory>
#include <thread>
#include <mutex>
#include <future>

namespace gwiz
{
	class VCFManager : public IVariantManager
	{
	public:
		typedef std::shared_ptr< VCFManager > SharedPtr;

		VCFManager(const std::vector< std::string >& vcfFilePaths, Region::SharedPtr regionPtr);
		~VCFManager();

		IVariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr) override;

		void asyncLoadVCFs();
		void waitForVCFsToLoadAndProcess();
		IVariantList::SharedPtr getCompleteVariantList() override;
	private:
		void processVCFs(); // a blocking call that waits for all vcfs to load and then combines them

		VariantList::SharedPtr m_variant_list_ptr;
		std::vector< std::string > m_vcf_paths;
		Region::SharedPtr m_region_ptr;
		std::atomic< bool > m_loaded_vcfs;
		std::mutex m_loaded_mutex;
		std::shared_ptr< std::thread > m_loading_thread_ptr;
	};
}

#endif //GWIZ_VCFMANAGER_H
