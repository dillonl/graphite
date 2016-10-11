#ifndef GRAPHITE_HTSLIB_ALIGNMENTREADER_HPP
#define GRAPHITE_HTSLIB_ALIGNMENTREADER_HPP

#include "IAlignmentReader.h"
#include "Sample.hpp"
#include "AlignmentReaderManager.hpp"

#include <mutex>

#include <htslib/sam.h>
#include <htslib/bgzf.h>

namespace graphite
{
	class HTSLibAlignmentReader : private IAlignmentReader, public std::enable_shared_from_this< HTSLibAlignmentReader >
	{
	public:
		typedef std::shared_ptr< HTSLibAlignmentReader > SharedPtr;

		HTSLibAlignmentReader(const std::string& filePath);
		HTSLibAlignmentReader(const std::string& filePath, AlignmentReaderManager< HTSLibAlignmentReader >* alignmentReaderManager);
		~HTSLibAlignmentReader();

		static std::vector< Sample::SharedPtr > GetBamReaderSamples(const std::string& bamPath);
		/* std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads = false); */
		std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads = false) override;
		void setAlignmentSequence(uint64_t filePosition, std::string& sequence);
		void open() override;
		void close() override;

		std::string getPath() { return m_file_path; }
		uint32_t getReaderID() { return m_id; }

		position getRegionLastPosition(Region::SharedPtr regionPtr);

	private:
		AlignmentReaderManager< HTSLibAlignmentReader >* m_alignment_reader_manager_ptr;

		std::mutex m_lock;
		std::string m_file_path;
		bool m_is_open;
		samFile* m_fp;
		bam_hdr_t* m_bam_header;
		bam1_t* m_bam_alignment;
		hts_idx_t* m_bam_index;
		static uint32_t s_hts_idx;
	};
}

#endif //GRAPHITE_HTSLIB_ALIGNMENTREADER_HPP
