#ifndef GWIZ_BAMALIGNMENTMANAGER_H
#define GWIZ_BAMALIGNMENTMANAGER_H

#include "IAlignmentManager.h"
#include "core/region/Region.h"

namespace gwiz
{
	class BamAlignmentManager : public IAlignmentManager
	{
	public:
		typedef std::shared_ptr< BamAlignmentManager > SharedPtr;
	    BamAlignmentManager(const std::string& bamPath, Region::SharedPtr regionPtr);
		~BamAlignmentManager();

		IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr);
		void asyncLoadAlignments();
		void waitForAlignmentsToLoad();
	private:
		std::string m_bam_path;
		Region::SharedPtr m_region_ptr;
	};
}

#endif //GWIZ_BAMALIGNMENTMANAGER_H
