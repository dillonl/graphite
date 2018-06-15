#ifndef GRAPHITE_BAMALIGNMENT_H
#define GRAPHITE_BAMALIGNMENT_H

#include "core2/util/Noncopyable.hpp"

#include "api/BamAlignment.h"

namespace graphite
{
	class BamAlignment : private Noncopyable
	{
	public:
		typedef std::shared_ptr< BamAlignment > SharedPtr;
		BamAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr);
		~BamAlignment();

		std::string getRead();
		std::string getReadName();
		Sample::SharedPtr getSamplePtr();

	private:
		Sample::SharedPtr m_sample_ptr;
		std::shared_ptr< BamTools::BamAlignment > m_bam_alignment_ptr;
	};
}

#endif //GRAPHITE_BAMALIGNMENT_H
