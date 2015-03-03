#include "BamAlignmentReaderManager.h"

namespace gwiz
{
	BamAlignmentReaderManager::BamAlignmentReaderManager(const std::string& bamPath) :
		m_bam_path(bamPath)
	{
	}

	BamAlignmentReaderManager::~BamAlignmentReaderManager()
	{
	}

	IAlignmentReader::SharedPtr BamAlignmentReaderManager::generateAlignmentReader()
	{
		return std::make_shared< gwiz::BamAlignmentReader >(this->m_bam_path);
	}
}
