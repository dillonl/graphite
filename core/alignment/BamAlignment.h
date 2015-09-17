#ifndef GRAPHITE_BAMALIGNMENT_H
#define GRAPHITE_BAMALIGNMENT_H

#include "core/region/Region.h"
#include "IAlignment.h"

#include "api/BamReader.h"
#include "api/BamAlignment.h"

typedef std::shared_ptr< BamTools::BamAlignment > BamAlignmentPtr;

namespace graphite
{
	class BamAlignment : public IAlignment
	{
	public:
		typedef std::shared_ptr< BamAlignment > SharedPtr;
		BamAlignment(BamAlignmentPtr bamAlignmentPtr);
        virtual ~BamAlignment();

		const char* getSequence() override { return m_sequence_string.c_str(); }
		const size_t getLength() override { return m_sequence_string.size(); };
		const position getPosition() override { return m_position; }
		const std::string getID() override { return m_id; }
		const bool isFirstMate() override {return m_first_mate;}
		const bool isMapped() override { return m_mapped; }
		const bool isReverseStrand() override { return m_reverse_strand; }
		const uint16_t getOriginalMapQuality() override { return m_original_map_quality; }
		/*
		const char* getSequence() override;
		const position getPosition() override;
		const size_t getLength() override;
		const std::string getID() override { return (this->m_bam_alignment_ptr->IsFirstMate()) ? this->m_bam_alignment_ptr->Name + "1" : this->m_bam_alignment_ptr->Name + "0"; }
		const bool isFirstMate() override {return this->m_bam_alignment_ptr->IsFirstMate();}
		const bool isMapped() override { return this->m_bam_alignment_ptr->IsMapped(); }
		const bool isReverseStrand() override { return this->m_bam_alignment_ptr->IsReverseStrand(); }
		const uint16_t getOriginalMapQuality() override { return this->m_bam_alignment_ptr->MapQuality; }
		*/

    private:
		char* m_sequence;
		size_t m_length;
		std::string m_sequence_string;
		position m_position;
		std::string m_id;
		bool m_first_mate;
		bool m_mapped;
		bool m_reverse_strand;
		uint16_t m_original_map_quality;
        /* BamAlignmentPtr m_bam_alignment_ptr; */
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
