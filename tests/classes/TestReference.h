#ifndef GWIZ_TESTS_CLASSES_TESTREFERENCE_H
#define GWIZ_TESTS_CLASSES_TESTREFERENCE_H

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
		    TestReference(std::string& sequence, Region::SharedPtr region) :
			    m_sequence(sequence), IReference(region)
			{
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

#endif //GWIZ_TESTS_CLASSES_TESTREFERENCE_H
