#ifndef GWIZ_BAMALIGNMENTREADER_H
#define GWIZ_BAMALIGNMENTREADER_H

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include <boost/noncopyable.hpp>

#include <memory>

#include "core/util/Types.h"
#include "IAlignmentReader.h"
#include "core/region/Region.h"

namespace gwiz
{
	class BamAlignmentReader : public IAlignmentReader
	{
	public:
		typedef std::shared_ptr< BamAlignmentReader > SharedPtr;

		static int32_t ReaderCount;

	    BamAlignmentReader(const std::string& bamPath);
		virtual ~BamAlignmentReader();

		virtual void init() override;
		virtual void releaseReader() override;
		virtual size_t getAverageReadLength() override { return m_average_bam_read_length; }
		virtual bool getNextAlignment(IAlignment::SharedPtr& alignment) override;

		void setRegion(Region::SharedPtr region) override;
		virtual size_t getReadCount() override { throw "getReadCount not implemented for BamAlignmentReader, use preloaded reader for count"; }

	private:
		void setAverageBamReadLength();
		size_t m_average_bam_read_length;
		std::string m_bam_path;
		BamTools::BamReader m_bam_reader;

	};
}

#endif //GWIZ_IALIGNMENTREADER_H
