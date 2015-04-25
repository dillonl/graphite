#include "BamAlignmentReader.h"
#include "BamAlignment.h"

namespace gwiz
{

	BamAlignmentReader::BamAlignmentReader(const std::string& bamPath) :
		m_bam_path(bamPath), m_average_bam_read_length(0)
	{
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
		releaseReader();
	}

	void BamAlignmentReader::init()
	{
		if (!m_bam_reader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
		setAverageBamReadLength();
	}

	void BamAlignmentReader::releaseReader()
	{
		m_bam_reader.Close();
	}

	void BamAlignmentReader::setAverageBamReadLength()
	{
		IAlignment::SharedPtr alignmentPtr;
		if (getNextAlignment(alignmentPtr))
		{
			m_average_bam_read_length = alignmentPtr->getLength();
			this->m_bam_reader.Rewind();
		}
	}

	bool BamAlignmentReader::getNextAlignment(IAlignment::SharedPtr& alignmentPtr)
	{
		BamAlignmentPtr bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		bool returnValue = this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr);
		alignmentPtr = std::make_shared< BamAlignment >(bamAlignmentPtr);

		return returnValue;
	}

	void BamAlignmentReader::setRegion(Region::SharedPtr region)
	{
		this->m_bam_reader.LocateIndex();
		int refID = this->m_bam_reader.GetReferenceID(region->getReferenceID());
		this->m_bam_reader.SetRegion(refID, region->getStartPosition(), refID, region->getEndPosition());
		this->m_region = region;
	}

}
