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
        // Currently a vector ref sequences is stored with each bamAlignment. A more optimized solution may be to create only a single instance. Maybe a static vector.
	    BamAlignment(BamTools::BamAlignment& bamAlignment, std::shared_ptr< Sample > samplePtr, BamTools::RefVector refVector) :
				m_position(bamAlignment.Position),
				m_first_mate(bamAlignment.IsFirstMate()),
				m_mapped(bamAlignment.IsMapped()),
				m_reverse_strand(bamAlignment.IsReverseStrand()),
				m_original_map_quality(bamAlignment.MapQuality),
				m_id(bamAlignment.Name + std::to_string(bamAlignment.IsFirstMate())),
				m_length(bamAlignment.QueryBases.size()),
				m_duplicate(bamAlignment.IsDuplicate()),
				m_sequence_counter(0),
                
                m_name(bamAlignment.Name),                      // Read name.
                m_alignment_flag(bamAlignment.AlignmentFlag),   // Alignment bit-flag.
                m_ref_ID(bamAlignment.RefID),              // Wrong value is returned... ID number for reference sequence.
                //m_cigar_data(bamAlignment.CigarData),           // CIGAR operations for this alignment.
                m_mateID(bamAlignment.Qualities),               // ID number for the reference sequence where alignment's mate was aligned.
                m_mate_position(bamAlignment.MatePosition),     // Position (0-based) where alignment's mate starts.
                m_template_length(bamAlignment.InsertSize),     // Ovserved template length
                m_qualities(bamAlignment.Qualities),            // FASTQ qualities (ASCII, not numerical values).
                m_ref_vector(refVector)                         // Original reference name of the alignment.
		{
			m_sample_ptr = samplePtr;
			m_sequence = new char[bamAlignment.QueryBases.size() + 1];
			memcpy(m_sequence, bamAlignment.QueryBases.c_str(), bamAlignment.QueryBases.size() + 1); // the +1 is for the '\0' char
		}
        // Get bam alignment information.
        // QNAME - BamTools::BamAlignment.Name???
        // FLAG - BamTools::BamAlignment.AlignmentFlag
        // RNAME - BamTools::BamAlignment.RefID???
        // POS - BamTools::BamAlignment.Position
        //          Needs to be recalculated from GSSWGraph
        // MAPQ - BamTools::BamAlignment.MapQuality
        // CIGAR - BamTools::BamAlignment.CigarData???
        //          Will need to recalculated (look at Dillon's code).
        // RNEXT - BamTools::BamAlignment.MateRefID
        // PNEXT - BamTools::BamAlignment.MatPosition
        // TLEN - BamTools::BamAlignment.InsertSize???
        // SEQ - BamTools::BamAlignment.QueryBases???
        // QUAL - BamTools::BamAlignment.Qualities

		virtual ~BamAlignment() { delete[] m_sequence; }

		const char* getSequence() override { return m_sequence; }
		const size_t getLength() override { return m_length; };
		const position getPosition() override { return m_position; }
		const std::string getID() override { return m_id; }
		const bool isFirstMate() override {return m_first_mate;}
		const bool isMapped() override { return m_mapped; }
		const bool isReverseStrand() override { return m_reverse_strand; }
		const bool isDuplicate() override { return m_duplicate; }
		const uint16_t getOriginalMapQuality() override { return m_original_map_quality; }

		const void incrementReferenceCount()
		{
			std::lock_guard< std::mutex > l(m_sequence_mutex);
			++m_sequence_counter;
		}
		const void setSequence(char* seq, uint32_t len) override
		{
			std::lock_guard< std::mutex > l(m_sequence_mutex);
			m_sequence = seq;
			m_length = len;
		}

		const void removeSequence() override
		{
			std::lock_guard< std::mutex > l(m_sequence_mutex);
			--m_sequence_counter;
			if (m_sequence_counter <= 0)
			{
				m_sequence_counter = 0;
				m_length = 0;
			}
		}

        // New functions added.
        const std::string getName () { return m_name; }
        const uint32_t getAlignmentFlag () { return m_alignment_flag; }
        const int32_t getRefID () { return m_ref_ID; } // Wrong value returned.
        //void setNewPosition () {};
        //const std::vector< BamTools::CigarOp > getOriginalCigarData () { return m_cigar_data; }
        //void calculateNewCigarData () {};
        //const std::vector getNewCigarData () { return m_new_cigar_data; }
        const std::string getMateID () { return m_mateID; }
        const int32_t getMatePosition () { return m_mate_position; }
        const int32_t getTemplateLength () { return m_template_length; }
        const std::string getFastqQualities () { return m_qualities; }

        std::string getAlignmentRegion (int32_t refID)
        {
            return m_ref_vector[refID].RefName;
        }

    private:
		std::mutex m_sequence_mutex;
		uint16_t m_sequence_counter;
		char* m_sequence;
		size_t m_length;
		position m_position;
		std::string m_id;
		bool m_first_mate;
		bool m_mapped;
		bool m_reverse_strand;
		bool m_duplicate;
		uint16_t m_original_map_quality;

        std::string m_name;
        uint32_t m_alignment_flag;
        int32_t m_ref_ID;
        std::string m_mateID;
        int32_t m_mate_position;
        int32_t m_template_length;
        std::string m_qualities;
        position m_new_position;                // May not need this.
        BamTools::RefVector m_ref_vector;       // BamTools::BamReader::RefVector
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
