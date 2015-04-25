#ifndef GWIZ_BAMALIGNMENTREADERPRELOAD_H
#define GWIZ_BAMALIGNMENTREADERPRELOAD_H

#include "IAlignmentReader.h"
#include "BamAlignment.h"

namespace gwiz
{
	class BamAlignmentReaderPreload : public IAlignmentReader
	{
	public:
		typedef std::shared_ptr< BamAlignmentReaderPreload > SharedPtr;

		BamAlignmentReaderPreload(std::shared_ptr< std::vector< BamAlignment::SharedPtr > > alignmentsPtr);
		virtual ~BamAlignmentReaderPreload();

		virtual void setRegion(Region::SharedPtr region) override;
		virtual size_t getAverageReadLength() override;
		virtual bool getNextAlignment(IAlignment::SharedPtr& alignment) override;
		virtual size_t getReadCount() override { return m_alignments_ptr->size(); }

	private:
		void setAverageBamReadLength();
		void setIndexClosestToPosition(position pos, size_t& index, bool subtractLength);
		std::shared_ptr< std::vector< BamAlignment::SharedPtr > > m_alignments_ptr;
		size_t m_start_index;
		size_t m_current_index;
		size_t m_end_index;
		size_t m_average_bam_read_length;

	};
}

#endif //GWIZ_BAMALIGNMENTREADERPRELOAD_H
