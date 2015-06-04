#include "VCFManager.h"
#include "Variant.h"

#include "core/parser/VCFParser.hpp"
#include "core/util/ThreadPool.hpp"

#include <functional>

namespace gwiz
{
	VCFManager::VCFManager(const std::vector< std::string >& vcfFilePaths, Region::SharedPtr regionPtr) :
		m_loaded_vcfs(false),
		m_vcf_paths(vcfFilePaths),
		m_region_ptr(regionPtr)
	{
	}

	VCFManager::~VCFManager()
	{
	}

	void VCFManager::asyncLoadVCFs()
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		if (this->m_loaded_vcfs) { return; }
		this->m_loading_thread_ptr = std::make_shared< std::thread >(&VCFManager::processVCFs, this);
	}

	// the intention of this function is to process the vcfs in a blocking way
	void VCFManager::processVCFs()
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		std::vector< std::thread > vcfLoadThreads;
		std::vector< std::shared_ptr< std::future< std::vector< IVariant::SharedPtr > > > > vcfFutureVariantListPtrs;
		for (const auto& vcfFilePath : this->m_vcf_paths)
		{
			auto vcfListPtr = std::make_shared< VCFList >(vcfFilePath, this->m_region_ptr);
			auto funct = std::bind(&VCFList::load, vcfListPtr);
			auto functFuture = ThreadPool::Instance()->enqueue(funct);
			vcfFutureVariantListPtrs.emplace_back(functFuture);
		}

		std::vector< IVariant::SharedPtr > variantPtrs;
		for (auto& vcfFuture : vcfFutureVariantListPtrs)
		{
			vcfFuture->wait();
			auto vcfVariantPtrs = vcfFuture->get();
			variantPtrs.insert(variantPtrs.end(), vcfVariantPtrs.begin(), vcfVariantPtrs.end());
		}
		vcfFutureVariantListPtrs.clear();
		this->m_variant_list_ptr = std::make_shared< VariantList >(variantPtrs);
		this->m_variant_list_ptr->sort();
		this->m_variant_list_ptr->normalizeOverlappingVariants();
		this->m_loaded_vcfs = true;
	}

	void VCFManager::waitForVCFsToLoadAndProcess()
	{
		this->m_loading_thread_ptr->join();
	}

	IVariantList::SharedPtr VCFManager::getVariantsInRegion(Region::SharedPtr regionPtr)
	{
		return this->m_variant_list_ptr->getVariantsInRegion(regionPtr);
	}

	IVariantList::SharedPtr VCFManager::getCompleteVariantList()
	{
		return this->m_variant_list_ptr;
	}
}
