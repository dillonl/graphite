#include "AlignmentList.h"
#include "SampleManager.hpp"
#include "SamtoolsAlignmentReader.h"

namespace graphite
{
	AlignmentList::AlignmentList(std::vector< IAlignment::SharedPtr > alignmentPtrs) :
		m_alignment_ptrs(alignmentPtrs),
		m_alignment_idx(0)
	{
	}

	AlignmentList::AlignmentList(std::vector< IAlignment::SharedPtr > alignmentPtrs, std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > nameAlignmentPtrMapPtr, Region::SharedPtr regionPtr) :
		m_alignment_ptrs(alignmentPtrs),
		m_alignment_idx(0),
		m_name_alignment_ptr_map_ptr(nameAlignmentPtrMapPtr)
	{
		uint32_t diff = 0;
		// auto startPos = (regionPtr->getStartPosition() - diff < 0) ? 0 : regionPtr->getStartPosition() - diff;
		auto startPos = alignmentPtrs[0]->getPosition() - 20;
		auto endPos = alignmentPtrs[alignmentPtrs.size() - 1]->getPosition() + 200;
		m_region_ptr = std::make_shared< Region >(regionPtr->getReferenceID(), startPos, endPos);
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

	void AlignmentList::loadAlignmentSequences()
	{
		for (auto& iter : SampleManager::Instance()->getSamplePtrs())
		{
			auto samplePtr = iter.second;
			auto readerPtr = std::make_shared< SamtoolsAlignmentReader >(samplePtr->getPath());
			readerPtr->loadAlignmentSequencesInRegion(this->m_region_ptr, m_name_alignment_ptr_map_ptr);
		}
		/*
		for (auto alignmentPtr : this->m_alignment_ptrs)
		{
			std::cout << alignmentPtr->getPosition() << " " << alignmentPtr->getSequence() << std::endl;
		}
		*/
	}

	void AlignmentList::unloadAlignmentSequences()
	{
		for (auto alignmentPtr : this->m_alignment_ptrs)
		{
			std::cout << "---removing---" << std::endl;
			alignmentPtr->removeSequence();
		}
	}
}
