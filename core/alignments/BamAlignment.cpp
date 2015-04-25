#include "BamAlignment.h"

namespace gwiz
{
	BamAlignment::BamAlignment(BamAlignmentPtr bamAlignmentPtr) :
		m_bam_alignment_ptr(bamAlignmentPtr)
	{
	}

	BamAlignment::~BamAlignment()
	{
	}

	const char* BamAlignment::getSequence()
	{
		return this->m_bam_alignment_ptr->QueryBases.c_str();
	}

	const position BamAlignment::getPosition()
	{
		return this->m_bam_alignment_ptr->Position + 1; // vcf and fasta are 1 based, bamtools is 0 based so we must add a 1
	}

	const size_t BamAlignment::getLength()
	{
		return this->m_bam_alignment_ptr->QueryBases.size();
	}
}
