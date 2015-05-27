#ifndef GWIZ_BAMALIGNMENTREADERMANAGER_H
#define GWIZ_BAMALIGNMENTREADERMANAGER_H

#include "IAlignmentReaderManager.h"
#include "BamAlignmentReader.h"

namespace gwiz
{
	class BamAlignmentReaderManager : public IAlignmentReaderManager
	{
	public:
		typedef std::shared_ptr< BamAlignmentReaderManager > SharedPtr;

		BamAlignmentReaderManager(const std::string& bamPath);
		~BamAlignmentReaderManager();

		virtual IAlignmentReader::SharedPtr generateAlignmentReader() override;

	protected:
		std::string m_bam_path;
	};
}

#endif //GWIZ_BAMALIGNMENTREADERMANAGER_H
