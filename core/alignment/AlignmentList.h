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
		~AlignmentList();

		size_t getCount() override;
		void sort() override;
		bool getNextAlignment(IAlignment::SharedPtr& alignmentPtr) override;
	private:
		std::vector< IAlignment::SharedPtr > m_alignment_ptrs;
		size_t m_alignment_idx;
	};
}

#endif //GRAPHITE_ALIGNMENTLIST_H
