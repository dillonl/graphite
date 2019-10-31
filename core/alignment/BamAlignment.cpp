#include "BamAlignment.h"

namespace graphite
{

	BamAlignment::BamAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr) :
		m_bam_alignment_ptr(bamAlignmentPtr),
		m_sample_ptr(samplePtr)
	{
	}

	BamAlignment::~BamAlignment()
	{
	}

	std::string BamAlignment::getRead()
	{
		return this->m_bam_alignment_ptr->QueryBases;
	}

	std::string BamAlignment::getReadName()
	{
		return this->m_bam_alignment_ptr->Name;
	}

	Sample::SharedPtr BamAlignment::getSamplePtr()
	{
		return m_sample_ptr;
	}

}
