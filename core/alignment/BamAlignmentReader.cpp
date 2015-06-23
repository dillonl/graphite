#include "BamAlignmentReader.h"
#include "BamAlignment.h"
#include "AlignmentList.h"

namespace gwiz
{
	BamAlignmentReader::BamAlignmentReader(const std::string& filePath) :
		m_file_path(filePath)
	{
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
	}

	std::vector< IAlignment::SharedPtr > BamAlignmentReader::loadAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		if (!this->m_bam_reader.Open(this->m_file_path))
		{
			throw "Unable to open bam file";
		}

		this->m_bam_reader.LocateIndex();
		int refID = this->m_bam_reader.GetReferenceID(regionPtr->getReferenceID());
		// add 1 to the start and end positions because this is 0 based
		this->m_bam_reader.SetRegion(refID, regionPtr->getStartPosition(), refID, regionPtr->getEndPosition());

		auto bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		size_t counter = 0;
		std::vector< IAlignment::SharedPtr > alignmentPtrs;
		uint32_t count = 0;
		while(this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr))
		{
			if (bamAlignmentPtr->RefID != refID) { break; }
			alignmentPtrs.push_back(std::make_shared< BamAlignment >(bamAlignmentPtr));
			bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		}
		return alignmentPtrs;
	}
}
