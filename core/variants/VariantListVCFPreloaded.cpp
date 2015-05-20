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
		uint32_t count = 1;
		// while (vcfFileReader->getNextCompoundVariant(variantPtr))
		while (vcfFileReader->getNextVariant(variantPtr))
		{
			// add all variants unless region is set
			if (this->m_region_ptr != nullptr)
			{
				bool sameReferenceID = (std::strcmp(this->m_region_ptr->getReferenceID().c_str(), variantPtr->getChrom().c_str()) == 0);
				if (sameReferenceID && (this->m_region_ptr->getEndPosition() < variantPtr->getPosition())) // dont continue; if the variants are past the end position
				{
					break;
				}
				if (!sameReferenceID || (sameReferenceID && (variantPtr->getPosition() < this->m_region_ptr->getStartPosition())))
				{
					continue;
				}
			}
			// std::cout << "variant: ";
			this->m_variants_ptr_list.push_back(variantPtr);
			// std::cout << count++ << std::endl;
		}
	}
}
