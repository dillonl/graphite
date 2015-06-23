#include "VCFManager.h"
#include "Variant.h"

#include "core/parser/VCFParser.hpp"
#include "core/util/ThreadPool.hpp"

#include <functional>

namespace gwiz
{
	VCFManager::VCFManager(const std::string& vcfPath, Region::SharedPtr regionPtr) :
		m_loaded_vcfs(false),
		m_region_ptr(regionPtr)
	{
		auto vcfFileReaderPtr = std::make_shared< VCFFileReader >(vcfPath);
		this->m_vcf_file_reader_ptrs.emplace_back(vcfFileReaderPtr);
	}

	VCFManager::VCFManager(const std::vector< std::string >& vcfFilePaths, Region::SharedPtr regionPtr) :
		m_loaded_vcfs(false),
		m_region_ptr(regionPtr)
	{
		for (const auto& vcfPath : vcfFilePaths)
		{
			auto vcfFileReaderPtr = std::make_shared< VCFFileReader >(vcfPath);
			this->m_vcf_file_reader_ptrs.emplace_back(vcfFileReaderPtr);
		}
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
		for (const auto& vcfFileReaderPtr : this->m_vcf_file_reader_ptrs)
		{
			auto funct = std::bind(&VCFFileReader::getVariantsInRegion, vcfFileReaderPtr, this->m_region_ptr);
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

	void VCFManager::printToVCF(std::ostream& out)
	{
		this->m_variant_list_ptr->printToVCF(out);
	}

	void VCFManager::releaseResources()
	{
	}
}
