#include "VariantListVCFPreloaded.h"

namespace gwiz
{
	VariantListVCFPreloaded::VariantListVCFPreloaded(const std::string& path, Region::SharedPtr region) :
		m_processed(false)
	{
	}

	VariantListVCFPreloaded::VariantListVCFPreloaded(const std::string& path) :
		m_processed(false)
	{
	}

	VariantListVCFPreloaded::~VariantListVCFPreloaded()
	{
	}

	void VariantListVCFPreloaded::rewind()
	{
	}

	void VariantListVCFPreloaded::process()
	{
	}

	bool VariantListVCFPreloaded::getNextVariant(Variant::SharedPtr& variant)
	{
		return false;
	}
	bool VariantListVCFPreloaded::getPreviousVariant(Variant::SharedPtr& variant)
	{
		return false;
	}

	bool VariantListVCFPreloaded::getVariant(Variant::SharedPtr& variant, const uint32_t index)
	{
		return false;
	}
}
