#ifndef GRAPHITE_ALIGNMENTLIST_H
#define GRAPHITE_ALIGNMENTLIST_H

#include "IAlignment.h"
#include "IAlignmentList.h"

namespace graphite
{
	class AlignmentList : public IAlignmentList
	{
	public:
		typedef std::shared_ptr< AlignmentList > SharedPtr;
		typedef std::vector< IAlignment::SharedPtr >::iterator AlignmentIter;
		AlignmentList(std::vector< IAlignment::SharedPtr > alignmentPtrs);
		AlignmentList(std::vector< IAlignment::SharedPtr > alignmentPtrs, std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > nameAlignmentPtrMapPtr, Region::SharedPtr regionPtr);
		~AlignmentList();

		size_t getCount() override;
		void sort() override;
		bool getNextAlignment(IAlignment::SharedPtr& alignmentPtr) override;
		void loadAlignmentSequences() override;
		void unloadAlignmentSequences() override;
	private:
		Region::SharedPtr m_region;
		std::vector< IAlignment::SharedPtr > m_alignment_ptrs;
		size_t m_alignment_idx;
		Region::SharedPtr m_region_ptr;
		std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > m_name_alignment_ptr_map_ptr;
	};
}

#endif //GRAPHITE_ALIGNMENTLIST_H
