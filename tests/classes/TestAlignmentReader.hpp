#ifndef GWIZ_TESTING_TESTALIGNMENTREADER_HPP
#define GWIZ_TESTING_TESTALIGNMENTREADER_HPP

#include "core/alignment/IAlignmentReader.h"
#include "TestAlignment.hpp"

#include <queue>

namespace gwiz
{
namespace testing
{
	class TestAlignmentReader : public IAlignmentReader
	{
	public:
		typedef std::shared_ptr< TestAlignmentReader > SharedPtr;
		TestAlignmentReader() : m_current_alignment_index(0)
		{
		}

		~TestAlignmentReader()
		{
		}

		void setRegion(const std::string& referenceID)
		{
			this->m_region = std::make_shared< Region >(referenceID + ":" + std::to_string(this->m_alignments.front()->getPosition()) + "-" + std::to_string(this->m_alignments.back()->getPosition()));
		}

		virtual size_t getAverageReadLength() override
		{
			size_t readLength = 0;
			if (!this->m_alignments.empty())
			{
				this->m_alignments.at(0)->getLength();
			}
			return readLength;
		}

		void addAlignment(const position pos, const std::string& seq)
		{
			auto alignment = std::make_shared< TestAlignment >(pos, seq);
			this->m_alignments.emplace_back(alignment);

		}

		bool getNextAlignment(IAlignment::SharedPtr& alignment)
		{
			if (this->m_current_alignment_index < this->m_alignments.size())
			{
				alignment = this->m_alignments.at(m_current_alignment_index);
				++m_current_alignment_index;
				return true;
			}
			else
			{
				return false;
			}
		}

		virtual size_t getReadCount() override { return m_alignments.size(); }

		std::vector< IAlignment::SharedPtr > m_alignments;
		size_t m_current_alignment_index;
	};
}
}

#endif // GWIZ_TESTING_TESTALIGNMENTREADER_HPP
