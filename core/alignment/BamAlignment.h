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
				m_length(bamAlignment.QueryBases.size())
		{
			m_sample_ptr = samplePtr;
			m_sequence = new char[bamAlignment.QueryBases.size() + 1];
			memcpy(m_sequence, bamAlignment.QueryBases.c_str(), bamAlignment.QueryBases.size() + 1); // the +1 is for the '\0' char
		}

		virtual ~BamAlignment() { delete[] m_sequence; }

		const char* getSequence() override { return m_sequence; }
		const size_t getLength() override { return m_length; };
		const position getPosition() override { return m_position; }
		const std::string getID() override { return m_id; }
		const bool isFirstMate() override {return m_first_mate;}
		const bool isMapped() override { return m_mapped; }
		const bool isReverseStrand() override { return m_reverse_strand; }
		const uint16_t getOriginalMapQuality() override { return m_original_map_quality; }

    private:
		char* m_sequence;
		size_t m_length;
		position m_position;
		std::string m_id;
		bool m_first_mate;
		bool m_mapped;
		bool m_reverse_strand;
		uint16_t m_original_map_quality;
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
