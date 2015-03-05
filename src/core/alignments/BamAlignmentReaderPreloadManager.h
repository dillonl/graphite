#ifndef GWIZ_BAMALIGNMENTREADERPRELOADMANAGER_H
#define GWIZ_BAMALIGNMENTREADERPRELOADMANAGER_H

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "IAlignmentReaderManager.h"
#include "BamAlignmentReaderPreload.h"
#include "BamAlignmentReader.h"
#include "BamAlignment.h"

namespace gwiz
{
	class BamAlignmentReaderPreloadManager : public IAlignmentReaderManager
	{
	public:
		typedef std::shared_ptr< BamAlignmentReaderPreloadManager > SharedPtr;
		BamAlignmentReaderPreloadManager(const std::string& bamPath);
		BamAlignmentReaderPreloadManager(const std::string& bamPath, Region::SharedPtr region);
		virtual ~BamAlignmentReaderPreloadManager();

		virtual IAlignmentReader::SharedPtr generateAlignmentReader() override;

	private:
		void processBam();

		std::shared_ptr< std::vector< BamAlignment::SharedPtr > > m_alignments_ptr;
		std::string m_bam_path;
		Region::SharedPtr m_region;
        BamTools::BamReader m_bam_reader;
	};
}

#endif //GWIZ_BAMALIGNMENTREADERPRELOADMANAGER_H
