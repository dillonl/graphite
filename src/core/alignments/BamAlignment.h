#ifndef GWIZ_BAMALIGNMENT_H
#define GWIZ_BAMALIGNMENT_H

#include "core/region/Region.h"
#include "IAlignment.h"

#include "api/BamReader.h"
#include "api/BamAlignment.h"

typedef std::shared_ptr< BamTools::BamAlignment > BamAlignmentPtr;

namespace gwiz
{
	class BamAlignment : public IAlignment
	{
	public:
		typedef std::shared_ptr< BamAlignment > SharedPtr;
		BamAlignment(BamAlignmentPtr bamAlignmentPtr);
        virtual ~BamAlignment();

		const char* getSequence() override;
		const position getPosition() override;
		const size_t getLength() override;
		const std::string getID() override {return this->m_bam_alignment_ptr->Name;}
		const bool isFirstMate() override {return this->m_bam_alignment_ptr->IsFirstMate();}
		const bool isMapped() override { return this->m_bam_alignment_ptr->IsMapped(); }
		const bool isReverseStrand() override { return this->m_bam_alignment_ptr->IsReverseStrand(); }

    private:
        BamAlignmentPtr m_bam_alignment_ptr;
	};
}

#endif //GWIZ_BAMALIGNMENT_H
