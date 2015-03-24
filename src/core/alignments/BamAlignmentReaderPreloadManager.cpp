#include "BamAlignmentReaderPreloadManager.h"

namespace gwiz
{
	BamAlignmentReaderPreloadManager::BamAlignmentReaderPreloadManager(const std::string& bamPath) :
		m_bam_path(bamPath), m_region(nullptr), m_alignments_ptr(std::make_shared< std::vector< BamAlignment::SharedPtr > >())
	{

		std::cout << this->m_alignments_ptr.get() << std::endl;
		processBam();
	}

	BamAlignmentReaderPreloadManager::BamAlignmentReaderPreloadManager(const std::string& bamPath, Region::SharedPtr region) :
		m_bam_path(bamPath), m_region(region), m_alignments_ptr(std::make_shared< std::vector< BamAlignment::SharedPtr > >())
	{
		processBam();
	}

	BamAlignmentReaderPreloadManager::~BamAlignmentReaderPreloadManager()
	{
	}

	IAlignmentReader::SharedPtr BamAlignmentReaderPreloadManager::generateAlignmentReader()
	{
		return std::make_shared< BamAlignmentReaderPreload >(this->m_alignments_ptr);
	}

	void BamAlignmentReaderPreloadManager::processBam()
	{
		if (!this->m_bam_reader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
		if (this->m_region != nullptr)
		{
			this->m_bam_reader.LocateIndex();
			int refID = this->m_bam_reader.GetReferenceID(this->m_region->getReferenceID());
			// add 1 to the start and end positions because this is 0 based
			this->m_bam_reader.SetRegion(refID, this->m_region->getStartPosition(), refID, this->m_region->getEndPosition());
		}
		BamAlignmentPtr bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		size_t counter = 0;
		while(this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr))
		{
			// if (counter++ > 1000) { break; }
			this->m_alignments_ptr->push_back(std::make_shared< BamAlignment >(bamAlignmentPtr));
			bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		}
	}
}
