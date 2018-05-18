#include "VCFManager.h"
#include "Variant.h"

#include "core/util/ThreadPool.hpp"

#include <functional>
#include <unordered_map>

namespace graphite
{
	VCFManager::VCFManager(const std::string& vcfPath, Region::SharedPtr regionPtr, IReference::SharedPtr referencePtr, uint32_t readLength) :
		m_loaded_vcfs(false),
		m_region_ptr(regionPtr),
		m_reference_ptr(referencePtr)
	{
		auto vcfFileReaderPtr = VCFFileReader::CreateVCFFileReader(vcfPath, referencePtr, readLength);
		this->m_vcf_file_reader_ptrs.emplace_back(vcfFileReaderPtr);
	}

	VCFManager::VCFManager(const std::vector< std::string >& vcfFilePaths, Region::SharedPtr regionPtr, IReference::SharedPtr referencePtr, uint32_t readLength) :
		m_loaded_vcfs(false),
		m_region_ptr(regionPtr),
		m_reference_ptr(referencePtr)
	{
		for (const auto& vcfPath : vcfFilePaths)
		{
			auto vcfFileReaderPtr = VCFFileReader::CreateVCFFileReader(vcfPath, referencePtr, readLength);
			this->m_vcf_file_reader_ptrs.emplace_back(vcfFileReaderPtr);
		}
	}

	VCFManager::~VCFManager()
	{
	}

	// the intention of this function is to process the vcfs in a blocking way
	void VCFManager::loadVariants()
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
			auto tmpVCFVariantListPtr = std::make_shared< VariantList >(vcfVariantPtrs, m_reference_ptr);
			tmpVCFVariantListPtr->processOverlappingAlleles();
			this->m_path_vcf_variant_list_ptrs_map.emplace(vcfFileReaderPtr, tmpVCFVariantListPtr);
			variantPtrs.insert(variantPtrs.end(), vcfVariantPtrs.begin(), vcfVariantPtrs.end());
		}

		vcfFutureVariantListPtrsMap.clear();
		this->m_variant_list_ptr = std::make_shared< VariantList >(variantPtrs, m_reference_ptr);
		this->m_variant_list_ptr->sort();
		this->m_variant_list_ptr->normalizeOverlappingVariants();
		this->m_loaded_vcfs = true;
	}

	std::unordered_map< VCFFileReader::SharedPtr, VariantList::SharedPtr > VCFManager::getVCFReadersAndVariantListsMap()
	{
		return this->m_path_vcf_variant_list_ptrs_map;
	}

	VariantList::SharedPtr VCFManager::getVariantsInRegion(Region::SharedPtr regionPtr)
	{
		return this->m_variant_list_ptr->getVariantsInRegion(regionPtr);
	}

	VariantList::SharedPtr VCFManager::getCompleteVariantList()
	{
		return this->m_variant_list_ptr;
	}
}
