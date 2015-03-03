#ifndef GWIZ_TESTING_TESTALIGNMENTREADER_HPP
#define GWIZ_TESTING_TESTALIGNMENTREADER_HPP

#include "core/alignments/IAlignmentReader.h"
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
		TestAlignmentReader()
		{
		}

		~TestAlignmentReader()
		{
		}

		void setRegion(Region::SharedPtr region) {}

		virtual size_t getAverageReadLength() override
		{
			size_t readLength = 0;
			if (!this->m_alignments.empty())
			{
				this->m_alignments.front()->getLength();
			}
			return readLength;
		}

		void addAlignment(const position pos, const std::string& seq)
		{
			auto alignment = std::make_shared< TestAlignment >(pos, seq);
			this->m_alignments.push(alignment);

		}

		bool getNextAlignment(IAlignment::SharedPtr& alignment)
		{
			if (!this->m_alignments.empty())
			{
				alignment = this->m_alignments.front();
				this->m_alignments.pop();
			}
			return !this->m_alignments.empty();
		}

		std::queue< IAlignment::SharedPtr > m_alignments;
	};
}
}

#endif // GWIZ_TESTING_TESTALIGNMENTREADER_HPP
