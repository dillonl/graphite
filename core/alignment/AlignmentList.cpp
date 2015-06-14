#include "AlignmentList.h"

namespace gwiz
{
	AlignmentList::AlignmentList(std::vector< IAlignment::SharedPtr > alignmentPtrs) :
		m_alignment_ptrs(alignmentPtrs),
		m_alignment_idx(0)
	{
	}

	AlignmentList::~AlignmentList()
	{
	}

	size_t AlignmentList::getCount()
	{
		return m_alignment_ptrs.size();
	}

	void AlignmentList::sort()
	{
	}

	bool AlignmentList::getNextAlignment(IAlignment::SharedPtr& alignmentPtr)
	{
		if (this->m_alignment_idx >= this->m_alignment_ptrs.size())
		{
			alignmentPtr = nullptr;
			return false;
		}
		else
		{
			alignmentPtr = this->m_alignment_ptrs[this->m_alignment_idx];
			this->m_alignment_idx += 1;
			return true;
		}
	}
}
