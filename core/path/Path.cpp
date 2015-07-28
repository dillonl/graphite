#include "Path.h"

namespace gwiz
{
	Path::Path()
	{
	}

	Path::~Path()
	{
	}

	std::vector< IAllele::SharedPtr > Path::getAllelePath()
	{
		return this->m_allele_ptrs;
	}

	uint32_t Path::getPathSWPercent()
	{
		return this->m_sw_percentage;
	}

	IAlignment::SharedPtr Path::getAlignment()
	{
		return this->m_alignment_ptr;
	}

	void Path::addAlleleToPath(IAllele::SharedPtr allelePtr)
	{
		this->m_allele_ptrs.emplace_back(allelePtr);
	}

	void Path::setPathSWPercent(uint32_t swPercent)
	{
		this->m_sw_percentage = swPercent;
	}

	void Path::setAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		this->m_alignment_ptr = alignmentPtr;
	}
}
