#ifndef GRAPHITE_HTSLIB_ALIGNMENT_MANAGER_H
#define GRAPHITE_HTSLIB_ALIGNMENT_MANAGER_H

#include "IAlignmentManager.h"
#include "AlignmentReaderManager.hpp"
#include "HTSLibAlignmentReader.h"
#include "core/variant/IVariantManager.h"
#include "Sample.hpp"

#include <thread>

namespace graphite
{
	class HTSLibAlignmentManager : public IAlignmentManager
	{
	public:
		typedef std::shared_ptr< HTSLibAlignmentManager > SharedPtr;
		HTSLibAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr, AlignmentReaderManager< HTSLibAlignmentReader >::SharedPtr alignmentReaderManagerPtr, bool excludeDuplicateReads = false);
		~HTSLibAlignmentManager();

        IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr);
		void releaseResources();
		void processMappingStatistics();
		void loadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);
		std::vector< std::shared_ptr< Sample > > getSamplePtrs() { return m_sample_ptrs; }

	private:
		void loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding);

		std::mutex m_loaded_mutex;
        bool m_exclude_duplicate_reads;
		std::vector< Sample::SharedPtr > m_sample_ptrs;
		std::mutex m_alignment_ptrs_lock;
		bool m_loaded;
		AlignmentReaderManager< HTSLibAlignmentReader >::SharedPtr m_alignment_reader_manager;
	};
}

#endif //GRAPHITE_HTSLIB_ALIGNMENT_MANAGER_H
