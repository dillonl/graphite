#include "BamAlignment.h"

namespace graphite
{
	BamAlignment::BamAlignment(BamAlignmentPtr bamAlignmentPtr) :
		m_position(bamAlignmentPtr->Position),
		m_sequence_string(bamAlignmentPtr->QueryBases),
		m_first_mate(bamAlignmentPtr->IsFirstMate()),
		m_mapped(bamAlignmentPtr->IsMapped()),
		m_reverse_strand(bamAlignmentPtr->IsReverseStrand()),
		m_original_map_quality(bamAlignmentPtr->MapQuality)
	{
		m_id = (bamAlignmentPtr->IsFirstMate()) ? bamAlignmentPtr->Name + "1" : bamAlignmentPtr->Name + "0";

	}

	BamAlignment::~BamAlignment()
	{
	}


	/*
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
	*/
}
