#ifndef GRAPHITE_BAMREADER_H
#define GRAPHITE_BAMREADER_H

#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"
#include "core/sample/Sample.h"

#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"

#include <memory>
#include <vector>
#include <string>
#include <unordered_set>

namespace graphite
{
    typedef BamTools::BamAlignment BamAlignment;
	class BamReader : private Noncopyable
	{
	public:
		typedef std::shared_ptr< BamReader > SharedPtr;
		BamReader(const std::string& filename);
		~BamReader();

		void fetchBamAlignmentPtrsInRegion(std::vector< std::shared_ptr< BamAlignment > >& bamAlignmentPtrs,  Region::SharedPtr regionPtr, bool unmappedOnly, bool includeDuplicateReads);

        std::unordered_set< Sample::SharedPtr > getSamplePtrs();
		uint32_t getReadLength();

	private:
		void initializeSamplePtrs();
		std::shared_ptr< BamTools::BamReader > m_bam_reader;
		std::unordered_set< Sample::SharedPtr > m_sample_ptrs;
		std::unordered_set< std::string > m_sample_names;
		std::string m_bam_path;
	};
}

#endif// GRAPHITE_BAMREADER_H
