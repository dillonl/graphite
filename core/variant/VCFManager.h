#ifndef GRAPHITE_VCFMANAGER_H
#define GRAPHITE_VCFMANAGER_H

#include "core/region/Region.h"
#include "VariantList.h"
#include "VCFFileReader.h"
#include "IVariantManager.h"
#include "IVariant.h"

#include <atomic>
#include <memory>
#include <thread>
#include <mutex>
#include <future>

namespace graphite
{
	class VCFManager : public IVariantManager
	{
	public:
		typedef std::shared_ptr< VCFManager > SharedPtr;

		VCFManager(const std::string& vcfPath, Region::SharedPtr regionPtr);
		VCFManager(const std::vector< std::string >& vcfFilePaths, Region::SharedPtr regionPtr);
		~VCFManager();

		IVariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr) override;

		void asyncLoadVCFs();
		void waitForVCFsToLoadAndProcess();
		IVariantList::SharedPtr getCompleteVariantList() override;
		void releaseResources() override;
		void printToVCF(std::ostream& out);
		std::unordered_map< std::string, VariantList::SharedPtr > getVCFPathsAndVariantListsMap();
	private:
		void processVCFs(); // a blocking call that waits for all vcfs to load and then combines them

		VariantList::SharedPtr m_variant_list_ptr;
		std::vector< VCFFileReader::SharedPtr > m_vcf_file_reader_ptrs;
		std::unordered_map< std::string, VariantList::SharedPtr > m_path_vcf_variant_list_ptrs_map; // these are the original vcf variant lists. They are used to reconstruct the original vcf file output
		Region::SharedPtr m_region_ptr;
		std::atomic< bool > m_loaded_vcfs;
		std::mutex m_loaded_mutex;
		std::shared_ptr< std::thread > m_loading_thread_ptr;
	};
}

#endif //GRAPHITE_VCFMANAGER_H
