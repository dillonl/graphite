#include "BamAlignmentReaderPreload.h"

namespace gwiz
{
	BamAlignmentReaderPreload::BamAlignmentReaderPreload(std::shared_ptr< std::vector< BamAlignment::SharedPtr > > alignmentsPtr) :
		m_alignments_ptr(alignmentsPtr),
		m_current_index(0),
		m_start_index(0),
		m_end_index(alignmentsPtr->size() - 1),
		m_average_bam_read_length(0)
	{
	}

	BamAlignmentReaderPreload::~BamAlignmentReaderPreload()
	{
	}

	size_t BamAlignmentReaderPreload::getAverageReadLength()
	{
		return 100;
		/*
		if (m_average_bam_read_length == 0)
		{
			setAverageBamReadLength();
		}
		return m_average_bam_read_length;
		*/
	}

	void BamAlignmentReaderPreload::setAverageBamReadLength()
	{
		auto oldIndex = this->m_current_index;
		IAlignment::SharedPtr alignmentPtr;
		if (getNextAlignment(alignmentPtr))
		{
			m_average_bam_read_length = alignmentPtr->getLength();
			this->m_current_index = oldIndex;
		}
	}

	void BamAlignmentReaderPreload::setRegion(Region::SharedPtr regionPtr)
	{
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();

		auto lowerBound = std::lower_bound(this->m_alignments_ptr->begin(), this->m_alignments_ptr->end(), nullptr, [startPosition](const BamAlignment::SharedPtr& alignmentPtr, const BamAlignment::SharedPtr& ignore) {
				return startPosition > alignmentPtr->getPosition();
			});
		auto upperBound = std::upper_bound(this->m_alignments_ptr->begin(), this->m_alignments_ptr->end(), nullptr, [endPosition](const BamAlignment::SharedPtr& ignore, const BamAlignment::SharedPtr& alignmentPtr) {
				return alignmentPtr->getPosition() > endPosition;
			});

		this->m_start_index = lowerBound - this->m_alignments_ptr->begin();
		this->m_end_index = upperBound - this->m_alignments_ptr->begin();
		this->m_current_index = this->m_start_index;
		this->m_region = regionPtr;
	}

	bool BamAlignmentReaderPreload::getNextAlignment(IAlignment::SharedPtr& alignmentPtr)
	{
		// std::cout << "getNextAlignment" << std::endl;
		if (this->m_current_index <= this->m_end_index)
		{
			// std::cout << "getNextAlignment true" << std::endl;
			alignmentPtr = this->m_alignments_ptr->at(this->m_current_index);
		    ++this->m_current_index;
			return true;
		}
		else
		{
			// std::cout << "getNextAlignment false" << std::endl;
			alignmentPtr = nullptr;
			return false;
		}
	}

	/*
	 * NOTE: Currently this only works with just one chromosome
	 */
	void BamAlignmentReaderPreload::setIndexClosestToPosition(position pos, size_t& index, bool subtractLength)
	{
		size_t startIndex = 0;
		size_t lastIndex = this->m_alignments_ptr->size() - 1;
		while (startIndex <= lastIndex)
		{
			size_t midIndex = (startIndex + lastIndex) / 2;
			auto midPosition = (subtractLength) ? this->m_alignments_ptr->at(midIndex)->getPosition() - this->m_alignments_ptr->at(midIndex)->getLength() : this->m_alignments_ptr->at(midIndex)->getPosition();
			if (pos > midPosition)
			{
				startIndex = midIndex + 1;
			}
			else if (pos < midPosition)
			{
				lastIndex = midIndex - 1;
			}
			else
			{
				index = midIndex;
				return;
			}
		}
		index = lastIndex;
	}
}
