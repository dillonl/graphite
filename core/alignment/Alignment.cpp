#include "Alignment.h"

#include "htslib/sam.h"

namespace graphite
{
	Alignment::Alignment(char* htslibSeq, uint32_t len, std::string& readName, bool forwardStrand, bool firstMate, uint16_t mapQuality, Sample::SharedPtr samplePtr) : m_len(len), m_read_name(readName), m_is_forward_strand(forwardStrand), m_is_first_mate(firstMate), m_map_quality(mapQuality), m_sample_ptr(samplePtr)
	{
		m_unique_read_name = (m_read_name) + std::to_string(firstMate);

		this->m_seq = (char*)malloc(len+1);
		int n;
		for(n=0; n < len; n++)
		{
			this->m_seq[n] = seq_nt16_str[bam_seqi(htslibSeq,n)];
		}
		this->m_seq[n] = 0;
	}

	Alignment::~Alignment()
	{
		free(this->m_seq);
	}
}
