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
        // Currently a refVector  is stored with each bamAlignment. A more optimized solution may be to create only a single instance. Maybe a static vector.
	    BamAlignment(BamTools::BamAlignment& bamAlignment, std::shared_ptr< Sample > samplePtr, BamTools::RefVector refVector) :
                m_name(bamAlignment.Name),                      // QNAME.   Read name.
				m_position(bamAlignment.Position),              // POS.     Start position.
                m_ref_id(bamAlignment.RefID),                   // Reference ID corresponding to the reference name in the header.
                m_ref_vector(refVector),                        // A vector of the original reference names corresponding to the reference IDs.
                m_alignment_flag(bamAlignment.AlignmentFlag),   // FLAG.    Alignment bit-flag.
				m_original_map_quality(bamAlignment.MapQuality),// MAPQ.    Original mapping quality.
                m_cigar_data(bamAlignment.CigarData),           // CIGAR.   CIGAR operations for this alignment.
                m_mate_ref_id(bamAlignment.MateRefID),          // Mate reference ID corresponding to the alignments mate.
                m_mate_position(bamAlignment.MatePosition),     // PNEXT. Position (0-based) where alignment's mate starts.
                m_template_length(bamAlignment.InsertSize),     // TLEN. Obvserved template length.
                m_qualities(bamAlignment.Qualities),            // FASTQ qualities (ASCII, not numerical values).
				m_first_mate(bamAlignment.IsFirstMate()),
				m_mapped(bamAlignment.IsMapped()),
				m_reverse_strand(bamAlignment.IsReverseStrand()),
				m_id(bamAlignment.Name + std::to_string(bamAlignment.IsFirstMate())),
				m_length(bamAlignment.QueryBases.size()),
				m_duplicate(bamAlignment.IsDuplicate()),
				m_sequence_counter(0)
                
		{
			m_sample_ptr = samplePtr;
			m_sequence = new char[bamAlignment.QueryBases.size() + 1];
			memcpy(m_sequence, bamAlignment.QueryBases.c_str(), bamAlignment.QueryBases.size() + 1); // the +1 is for the '\0' char
		}

		virtual ~BamAlignment() { delete[] m_sequence; }

        const std::string getName () { return m_name; }
        const uint32_t getAlignmentFlag () { return m_alignment_flag; }
        std::string getReferenceName () { return m_ref_vector[m_ref_id].RefName; }  // Return the reference name using the reference ID.
		const position getPosition() override { return m_position; }
		const uint16_t getOriginalMapQuality() override { return m_original_map_quality; }

        const std::string getCigarString ()         // Prints out the original CIGAR string as a string.
        {
            std::string cigarString;
            //for (int i = 0; i < m_cigar_data.size(); ++i)
            for (auto cigarOp : m_cigar_data)
            {
                cigarString += std::to_string(cigarOp.Length);
                cigarString += cigarOp.Type;
            }
            
            if (cigarString.empty())
                return "*";
            else
                return cigarString;
        }

        std::string getMateReferenceName ()         // Return the reference name of the alignment's mate.
        {
            if (m_ref_id == m_mate_ref_id)
                return "=";
            else
                return m_ref_vector[m_mate_ref_id].RefName;
        }  

        const int32_t getMatePosition () { return m_mate_position; }        // Return the mate's position.
        const int32_t getTemplateLength () { return m_template_length; }    // Return the template length.
		const char* getSequence() override { return m_sequence; }
        const std::string getFastqQualities () { return m_qualities; }      // Return base quality scores.

		const size_t getLength() override { return m_length; };
		const std::string getID() override { return m_id; }
		const bool isFirstMate() override {return m_first_mate;}
		const bool isMapped() override { return m_mapped; }
		const bool isReverseStrand() override { return m_reverse_strand; }
		const bool isDuplicate() override { return m_duplicate; }

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

        //void setNewPosition () {};
        //const std::vector< BamTools::CigarOp > getOriginalCigarData () { return m_cigar_data; }
        //void calculateNewCigarData () {};
        //const std::vector getNewCigarData () { return m_new_cigar_data; }

    private:
		std::mutex m_sequence_mutex;
        std::string m_name;                 // QNAME. Sequence name.
        uint32_t m_alignment_flag;          // FLAG. Algnment flags.
        int32_t m_ref_id;                   // Reference ID corresponding to the reference name in the header.
        BamTools::RefVector m_ref_vector;   // A vector of reference names corresponding to reference IDs.
		position m_position;                // POS. Alignment starting position.
		uint16_t m_original_map_quality;    // MAPQ. Original mapping quality.
        std::vector< BamTools::CigarOp > m_cigar_data;  // CIGAR. Original CIGAR string for alignment.
        int32_t m_mate_ref_id;              // Mate reference ID for the alignments mate.
        int32_t m_mate_position;            // PNEXT. Position (0-based) where alignment's mate starts.
        int32_t m_template_length;          // TLEN. Obvserved template length.
		char* m_sequence;
        std::string m_qualities;            // FASTQ qualities (ASCII, not numerical values).
		uint16_t m_sequence_counter;
		size_t m_length;
		std::string m_id;
		bool m_first_mate;
		bool m_mapped;
		bool m_reverse_strand;
		bool m_duplicate;

        position m_new_position;                // May not need this.
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
