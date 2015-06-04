#include "VCFList.h"

namespace gwiz
{
	VCFList::VCFList(const std::string& filePath, Region::SharedPtr regionPtr)
	{
	}

	VCFList::~VCFList()
	{
	}

	bool VCFList::getNextVariant(IVariant::SharedPtr& variantPtr)
	{
		// check for region in here
		return false;
	}

	bool VCFList::hasVariants()
	{
		return false;
	}

	/*
	//things to be removed possibly
	size_t VCFList::getCount()
	{
		return 0;
	}

	// void VCFList::sort()
	// {
	// }

	void VCFList::addVariants(IVariantList::SharedPtr variantListPtr)
	{
	}
	*/

	std::vector< IVariant::SharedPtr > VCFList::load()
	{
		// skip to the start of the region here
		std::vector< IVariant::SharedPtr > variantPtrs;
		IVariant::SharedPtr variantPtr;
		while (getNextVariant(variantPtr))
		{
			variantPtrs.emplace_back(variantPtr);
		}
		return variantPtrs;
	}

	void VCFList::printToVCF(std::ostream& out)
	{
	}
}
