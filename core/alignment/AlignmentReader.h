#pragma once

#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"
#include "core/sample/Sample.h"
#include "Alignment.h"

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"

#include <string>
#include <unordered_map>

namespace graphite
{
	class AlignmentReader : private Noncopyable
	{
	public:
		typedef std::shared_ptr< AlignmentReader > SharedPtr;
	    AlignmentReader(const std::string& filename, const std::string& refPath);
	    ~AlignmentReader();

		void overwriteSample(Sample::SharedPtr samplePtr);
        bool shouldOverwriteSample() { return m_overwrite_sample; }
        void fetchAlignmentPtrsInRegion(std::vector< std::shared_ptr< Alignment > >& alignmentPtrs, Region::SharedPtr regionPtr, bool unmappedOnly, bool includeDuplicateReads, int32_t mappingQuality);
		std::unordered_map< std::string, Sample::SharedPtr > getSamplePtrs() { return this->m_sample_ptrs; }
		uint32_t getReadLength() { return this->m_read_length; }

	private:
		void init();
		void setReadLength();

		htsFile* m_in;
		bam_hdr_t* m_header;
		hts_idx_t* m_idx;
		uint32_t m_read_length;
		bool m_overwrite_sample;
		std::string m_path;
		std::unordered_map< std::string, Sample::SharedPtr > m_sample_ptrs;
		std::vector< std::string > m_available_regions;
    };
}
