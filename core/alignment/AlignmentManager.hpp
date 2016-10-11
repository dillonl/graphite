#ifndef GRAPHITE_ALIGNMENT_MANAGER_H
#define GRAPHITE_ALIGNMENT_MANAGER_H

#include "AlignmentReaderManager.hpp"
#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"

namespace graphite
{
	template < class AlignmentReaderType >
	class AlignmentManager : private Noncopyable
	{
	public:
		typedef std::shared_ptr< AlignmentManager< AlignmentReaderType > > SharedPtr;
		AlignmentManager(std::vector< std::string >& fileNames, Region::SharedPtr regionPtr) :
			m_alignment_reader_manager_ptr(std::shared_ptr< AlignmentReaderType >(fileNames, 20))
		{
		}

		~AlignmentManager()
		{
		}

		void loadAlignments()
		{

		}

	private:
		std::shared_ptr< AlignmentReaderManager< AlignmentReaderType > > m_alignment_reader_manager_ptr;
	};
}

#endif //GRAPHITE_ALIGNMENT_MANAGER_H
