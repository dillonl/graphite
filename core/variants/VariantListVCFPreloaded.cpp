#include "VCFFileReader.h"
#include "VariantListVCFPreloaded.h"

namespace gwiz
{
	VariantListVCFPreloaded::VariantListVCFPreloaded(const std::string& path, Region::SharedPtr region) :
		m_processed(false), m_path(path), m_region_ptr(region)
	{
	}

	VariantListVCFPreloaded::VariantListVCFPreloaded(const std::string& path) :
		m_processed(false), m_path(path), m_region_ptr(nullptr)
	{
	}

	VariantListVCFPreloaded::~VariantListVCFPreloaded()
	{
	}

	void VariantListVCFPreloaded::loadVariantsFromFile()
	{
		// std::cout << "loading..." << std::endl;
		auto vcfFileReader = std::make_shared< VCFFileReader >(this->m_path);
		Variant::SharedPtr variantPtr;
		while (vcfFileReader->getNextVariant(variantPtr))
		{
			// std::cout << "variant position: " << variantPtr->getPosition() << " " << this->m_region_ptr->getEndPosition() << std::endl;
			// add all variants unless region is set
			if (this->m_region_ptr != nullptr)
			{
				// if the variant isn't the same referenceID or if the variant position is less than the start position then continue;
				if (std::strcmp(this->m_region_ptr->getReferenceID().c_str(), variantPtr->getChrom().c_str()) != 0 || (this->m_region_ptr->getStartPosition() > variantPtr->getPosition() && this->m_region_ptr->getStartPosition() > 0 && this->m_region_ptr->getEndPosition() > 0))
				{
					continue;
				}
				// dont continue; if the variants are past the end position
				if (variantPtr->getPosition() > this->m_region_ptr->getEndPosition())
				{
					break;
				}
			}
			this->m_variants_ptr_list.push_back(variantPtr);
		}
	}
}
