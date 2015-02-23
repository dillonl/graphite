#ifndef GWIZ_BAMALIGNMENTREADER_H
#define GWIZ_BAMALIGNMENTREADER_H

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include <boost/noncopyable.hpp>

#include <memory>

#include "core/utils/Types.h"
#include "IAlignmentReader.h"

namespace gwiz
{
	class BamAlignmentReader : public IAlignmentReader
	{
	public:
		typedef std::shared_ptr< BamAlignmentReader > SharedPtr;

	    BamAlignmentReader(const std::string& bamPath);
		virtual ~BamAlignmentReader();

		virtual bool getNextAlignment(IAlignment::SharedPtr alignment) override;

	private:
		std::string m_bam_path;
		BamTools::BamReader m_bam_reader;

	};
}

#endif //GWIZ_IALIGNMENTREADER_H
