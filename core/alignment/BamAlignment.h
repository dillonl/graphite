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

	    BamAlignment(std::shared_ptr< BamTools::BamAlignment > bamtoolsAlignmentPtr, std::shared_ptr< Sample > samplePtr) :
		            m_bamtools_alignment_ptr(bamtoolsAlignmentPtr),
					m_id(bamtoolsAlignmentPtr->Name + std::to_string(bamtoolsAlignmentPtr->IsFirstMate()))
	    {
			m_sample_ptr = samplePtr;
	    }

		virtual ~BamAlignment() { }

		const char* getSequence() override { return m_bamtools_alignment_ptr->QueryBases.c_str(); }
		const size_t getLength() override { return  m_bamtools_alignment_ptr->QueryBases.size(); };
		const position getPosition() override { return m_bamtools_alignment_ptr->Position; }
		const std::string getID() override { return m_id; }
		const bool isFirstMate() override {return m_bamtools_alignment_ptr->IsFirstMate();}
		const bool isMapped() override { return m_bamtools_alignment_ptr->IsMapped(); }
		const bool isReverseStrand() override { return m_bamtools_alignment_ptr->IsReverseStrand(); }
		const bool isDuplicate() override { return m_bamtools_alignment_ptr->IsDuplicate(); }
		const uint16_t getOriginalMapQuality() override { return m_bamtools_alignment_ptr->MapQuality; }
		const std::shared_ptr< BamTools::BamAlignment > getBamToolsAlignmentPtr() override { return m_bamtools_alignment_ptr; }

    private:
		std::string m_id;
		std::shared_ptr< BamTools::BamAlignment > m_bamtools_alignment_ptr;
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
