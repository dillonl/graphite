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
		/*
	    BamAlignment(BamAlignmentPtr bamAlignmentPtr, std::shared_ptr< Sample > samplePtr) :
		            m_position(bamAlignmentPtr->Position),
					m_first_mate(bamAlignmentPtr->IsFirstMate()),
					m_mapped(bamAlignmentPtr->IsMapped()),
					m_reverse_strand(bamAlignmentPtr->IsReverseStrand()),
					m_original_map_quality(bamAlignmentPtr->MapQuality),
					m_id(bamAlignmentPtr->Name + std::to_string(bamAlignmentPtr->IsFirstMate())),
					m_length(bamAlignmentPtr->QueryBases.size())
		{
			m_sample_ptr = samplePtr;
			m_sequence = new char[bamAlignmentPtr->QueryBases.size() + 1];
			memcpy(m_sequence, bamAlignmentPtr->QueryBases.c_str(), bamAlignmentPtr->QueryBases.size() + 1); // the +1 is for the '\0' char
		}
		*/

	    BamAlignment(BamTools::BamAlignment& bamAlignment, std::shared_ptr< Sample > samplePtr) :
				m_position(bamAlignment.Position),
				m_first_mate(bamAlignment.IsFirstMate()),
				m_mapped(bamAlignment.IsMapped()),
				m_reverse_strand(bamAlignment.IsReverseStrand()),
				m_original_map_quality(bamAlignment.MapQuality),
				m_id(bamAlignment.Name + std::to_string(bamAlignment.IsFirstMate())),
				m_length(bamAlignment.QueryBases.size()),
				m_duplicate(bamAlignment.IsDuplicate()),
				m_sequence_counter(0)
		{
			m_sample_ptr = samplePtr;
			m_sequence = new char[bamAlignment.QueryBases.size() + 1];
			memcpy(m_sequence, bamAlignment.QueryBases.c_str(), bamAlignment.QueryBases.size() + 1); // the +1 is for the '\0' char
			m_tlen = bamAlignment.InsertSize;
		}

				/*
	    BamAlignment(position pos, bool firstMate, bool isMapped, bool isReverseStrand, bool duplicate, uint16_t mapQuality, char* id, std::shared_ptr< Sample > samplePtr) :
				m_position(pos),
				m_first_mate(firstMate),
				m_mapped(isMapped),
				m_reverse_strand(isReverseStrand),
				m_original_map_quality(mapQuality),
				m_id(std::string(id) + std::to_string(firstMate)),
				m_duplicate(duplicate),
				m_sequence_counter(0),
				m_length(0),
				m_sequence(nullptr)
		{
			m_sample_ptr = samplePtr;
		}
				*/
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

		int32_t getTLen() override
		{
			return m_tlen;
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
		int32_t m_tlen;
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
