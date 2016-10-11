#ifndef GRAPHITE_BAMALIGNMENTREADER_HPP
#define GRAPHITE_BAMALIGNMENTREADER_HPP

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "core/region/Region.h"
#include "IAlignment.h"
#include "IAlignmentList.h"
#include "IAlignmentReader.h"
#include "Sample.hpp"


namespace graphite
{
	class BamAlignmentReader : public IAlignmentReader
	{
	public:
        typedef std::shared_ptr< BamAlignmentReader > SharedPtr;
        BamAlignmentReader(const std::string& bamPath);
		~BamAlignmentReader();

		void open() override;
		void close() override;

        std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads = false) override;

		static std::vector< Sample::SharedPtr > GetBamReaderSamples(const std::string& bamPath);
		static position GetLastPositionInBam(const std::string& bamPath, Region::SharedPtr regionPtr);

	private:
        BamTools::BamReader m_bam_reader;
		std::string m_bam_path;
		bool m_is_open;
		std::mutex m_lock;
	};
}

#endif //GRAPHITE_BAMALIGNMENTREADER_HPP
