#ifndef GWIZ_TESTS_CLASSES_TESTREFERENCE_HPP
#define GWIZ_TESTS_CLASSES_TESTREFERENCE_HPP

#include "core/reference/IReference.h"

namespace gwiz
{
	namespace testing
	{

		/*
		 * This class is only for testing.
		 */

		class TestReference : public IReference
		{
		public:
			typedef std::shared_ptr< TestReference > SharedPtr;
		    TestReference(const std::string& sequence, const std::string& regionString) :
			    m_sequence(sequence)
			{
				m_region = std::make_shared< Region >(regionString);
			}

			~TestReference()
			{
			}

			const char* getSequence() override { return m_sequence.c_str(); }
			size_t getSequenceSize() override { return m_sequence.size(); }

		private:
			std::string m_sequence;

		};

	}
}

#endif //GWIZ_TESTS_CLASSES_TESTREFERENCE_HPP
