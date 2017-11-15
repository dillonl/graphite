#ifndef GRAPHITE_BAMALIGNMENTREADER_HPP
#define GRAPHITE_BAMALIGNMENTREADER_HPP

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "core/region/Region.h"

#include "AlignmentReaderManager.hpp"
#include "IAlignment.h"
#include "IAlignmentList.h"
#include "IAlignmentReader.h"


namespace graphite
{
	class BamAlignmentReader : public IAlignmentReader, public std::enable_shared_from_this< BamAlignmentReader >
	{
	public:
        typedef std::shared_ptr< BamAlignmentReader > SharedPtr;
        BamAlignmentReader(const std::string& bamPath);
		BamAlignmentReader(const std::string& bamPath, AlignmentReaderManager< BamAlignmentReader >* alignmentReaderManager);
		~BamAlignmentReader();

		void open() override;
		void close() override;

        std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, SampleManager::SharedPtr sampleManagerPtr, bool includeDuplicateReads = false) override;

        static std::vector< Sample::SharedPtr > GetBamReaderSamples(const std::string& bamPath);
		static position GetLastPositionInBam(const std::string& bamPath, Region::SharedPtr regionPtr);
		static uint32_t GetReadLength(const std::string& bamPath);

		std::string getPath() { return m_bam_path; }
        uint32_t getReaderID() { return m_id; }

	private:
        AlignmentReaderManager< BamAlignmentReader >* m_alignment_reader_manager_ptr;

		std::shared_ptr< BamTools::BamReader > m_bam_reader;
		std::string m_bam_path;
		bool m_is_open;
		std::mutex m_lock;
	};
}

#endif //GRAPHITE_BAMALIGNMENTREADER_HPP
