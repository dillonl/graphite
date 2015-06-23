#include "BamAlignmentManager.h"

namespace gwiz
{
	BamAlignmentManager::BamAlignmentManager(const std::string& bamPath, Region::SharedPtr regionPtr) :
		m_bam_path(bamPath),
		m_region_ptr(regionPtr)
	{
	}

	BamAlignmentManager::~BamAlignmentManager()
	{
	}

	IAlignmentList::SharedPtr BamAlignmentManager::getAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		// check if outside of the m_region_ptr and if so then throw an error
		return nullptr;
	}

	void BamAlignmentManager::asyncLoadAlignments()
	{
		std::lock_guard< std::mutex > lock(this->m_loaded_mutex);
		if (this->m_loaded) { return; }
	}

	void BamAlignmentManager::waitForAlignmentsToLoad()
	{
	}

	void BamAlignmentManager::releaseResources()
	{
	}
}
