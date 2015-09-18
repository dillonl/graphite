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
	    BamAlignment(BamAlignmentPtr bamAlignmentPtr) :
		            m_position(bamAlignmentPtr->Position),
					m_sequence_string(bamAlignmentPtr->QueryBases),
					m_first_mate(bamAlignmentPtr->IsFirstMate()),
					m_mapped(bamAlignmentPtr->IsMapped()),
					m_reverse_strand(bamAlignmentPtr->IsReverseStrand()),
					m_original_map_quality(bamAlignmentPtr->MapQuality),
					m_id(bamAlignmentPtr->Name + std::to_string(bamAlignmentPtr->IsFirstMate()))
		{
		}
		virtual ~BamAlignment() {}

		const char* getSequence() override { return m_sequence_string.c_str(); }
		const size_t getLength() override { return m_sequence_string.size(); };
		const position getPosition() override { return m_position; }
		const std::string getID() override { return m_id; }
		const bool isFirstMate() override {return m_first_mate;}
		const bool isMapped() override { return m_mapped; }
		const bool isReverseStrand() override { return m_reverse_strand; }
		const uint16_t getOriginalMapQuality() override { return m_original_map_quality; }

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
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
