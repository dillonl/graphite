#ifndef GWIZ_BAMALIGNMENTREADER_HPP
#define GWIZ_BAMALIGNMENTREADER_HPP

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "core/region/Region.h"
#include "IAlignment.h"
#include "IAlignmentList.h"
#include "IAlignmentReader.h"

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class BamAlignmentReader : private IAlignmentReader
	{
	public:
        typedef std::shared_ptr< BamAlignmentReader > SharedPtr;
        BamAlignmentReader(const std::string& filePath);
		~BamAlignmentReader();

        std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr) override;

	private:
		std::string m_file_path;
        BamTools::BamReader m_bam_reader;
	};
}

#endif //GWIZ_BAMALIGNMENTREADER_HPP
