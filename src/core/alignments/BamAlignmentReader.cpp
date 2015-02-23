#include "BamAlignmentReader.h"
#include "BamAlignment.h"

namespace gwiz
{

	BamAlignmentReader::BamAlignmentReader(const std::string& bamPath) :
		m_bam_path(bamPath)
	{
		if (!m_bam_reader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
		m_bam_reader.Close();
	}

	bool BamAlignmentReader::getNextAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		auto bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		bool returnValue = this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr);
		alignmentPtr = std::make_shared< BamAlignment >(bamAlignmentPtr);
		return returnValue;
	}

}
