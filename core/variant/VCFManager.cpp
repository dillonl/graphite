#include "VCFManager.h"
#include "Variant.h"

#include "core/parser/VCFParser.hpp"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <unordered_map>

namespace graphite
{
	VCFManager::VCFManager(const std::string& vcfPath, Region::SharedPtr regionPtr) :
		m_loaded_vcfs(false),
		m_region_ptr(regionPtr)
	{
		auto vcfFileReaderPtr = VCFFileReader::CreateVCFFileReader(vcfPath);
		this->m_vcf_file_reader_ptrs.emplace_back(vcfFileReaderPtr);
	}

	VCFManager::VCFManager(const std::vector< std::string >& vcfFilePaths, Region::SharedPtr regionPtr) :
		m_loaded_vcfs(false),
		m_region_ptr(regionPtr)
	{
		for (const auto& vcfPath : vcfFilePaths)
		{
			auto vcfFileReaderPtr = VCFFileReader::CreateVCFFileReader(vcfPath);
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
		std::unordered_map< std::string, std::shared_ptr< std::future< std::vector< IVariant::SharedPtr > > > > vcfFutureVariantListPtrsMap;
		for (const auto& vcfFileReaderPtr : this->m_vcf_file_reader_ptrs)
		{
			auto funct = std::bind(&VCFFileReader::getVariantsInRegion, vcfFileReaderPtr, this->m_region_ptr);
			auto functFuture = ThreadPool::Instance()->enqueue(funct);
			vcfFutureVariantListPtrsMap.emplace(vcfFileReaderPtr->getFilePath(), functFuture);
		}

		std::vector< IVariant::SharedPtr > variantPtrs;
		for (const auto& vcfFileReaderPtr : this->m_vcf_file_reader_ptrs)
		{
			std::string filePath = vcfFileReaderPtr->getFilePath();
			auto  vcfFuturePtr = vcfFutureVariantListPtrsMap[filePath];
			vcfFuturePtr->wait();
			auto vcfVariantPtrs = vcfFuturePtr->get();
			auto tmpVCFVariantListPtr = std::make_shared< VariantList >(vcfVariantPtrs);
			tmpVCFVariantListPtr->processOverlappingAlleles();
			this->m_path_vcf_variant_list_ptrs_map.emplace(filePath, tmpVCFVariantListPtr);
			variantPtrs.insert(variantPtrs.end(), vcfVariantPtrs.begin(), vcfVariantPtrs.end());
		}

		vcfFutureVariantListPtrsMap.clear();
		this->m_variant_list_ptr = std::make_shared< VariantList >(variantPtrs);
		this->m_variant_list_ptr->sort();
		this->m_variant_list_ptr->normalizeOverlappingVariants();
		// for (auto iter : this->m_path_vcf_variant_list_ptrs_map)
		// {
		// 	iter.second->processOverlappingAlleles();
		// }
		this->m_loaded_vcfs = true;
	}

	std::unordered_map< std::string, VariantList::SharedPtr > VCFManager::getVCFPathsAndVariantListsMap()
	{
		return this->m_path_vcf_variant_list_ptrs_map;
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
