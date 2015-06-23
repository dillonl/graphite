#ifndef GWIZ_BAMALIGNMENTMANAGER_H
#define GWIZ_BAMALIGNMENTMANAGER_H

#include "IAlignmentManager.h"
#include "core/region/Region.h"

#include <mutex>

namespace gwiz
{
	class BamAlignmentManager : public IAlignmentManager
	{
	public:
		typedef std::shared_ptr< BamAlignmentManager > SharedPtr;
	    BamAlignmentManager(const std::string& bamPath, Region::SharedPtr regionPtr);
		~BamAlignmentManager();

		IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr) override;
		void asyncLoadAlignments();
		void waitForAlignmentsToLoad();
		void releaseResources() override;
	private:

		std::mutex m_loaded_mutex;
		bool m_loaded;
		std::string m_bam_path;
		Region::SharedPtr m_region_ptr;
	};
}

#endif //GWIZ_BAMALIGNMENTMANAGER_H
