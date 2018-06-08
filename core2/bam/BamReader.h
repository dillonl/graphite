#ifndef GRAPHITE_BAMREADER_H
#define GRAPHITE_BAMREADER_H

#include "core2/util/Noncopyable.hpp"
#include "core2/region/Region.h"
#include "core2/sample/Sample.h"

#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamAux.h"

#include <memory>
#include <vector>
#include <string>

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

        std::vector< Sample::SharedPtr > getSamplePtrs();
		uint32_t getReadLength();

	private:
		std::shared_ptr< BamTools::BamReader > m_bam_reader;
		std::string m_bam_path;
	};
}

#endif// GRAPHITE_BAMREADER_H
