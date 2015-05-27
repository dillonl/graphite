#ifndef GWIZ_TESTS_SEQUENCEMANAGERTESTS_HPP
#define GWIZ_TESTS_SEQUENCEMANAGERTESTS_HPP

#include "core/sequence/SequenceManager.h"

class SequenceManagerTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        gwiz::SequenceManager::Instance()->clearSequences();
    }
};

// Tests factorial of negative numbers.
TEST_F(SequenceManagerTest, GetSequence)
{
	const char* sequence = "ATGC";
	auto sequencePtr = gwiz::SequenceManager::Instance()->getSequence(sequence);
	ASSERT_STREQ(sequence, sequencePtr->getSequence());
}

TEST_F(SequenceManagerTest, CheckSameSequenceSamePtr)
{
	const char* sequence = "ATGC";
	auto sequencePtrFirst = gwiz::SequenceManager::Instance()->getSequence(sequence);
	auto sequencePtrSecond = gwiz::SequenceManager::Instance()->getSequence(sequence);
	ASSERT_EQ(sequencePtrFirst.get(), sequencePtrSecond.get());
}

TEST_F(SequenceManagerTest, CheckDifferentSequenceDifferentPtr)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
	auto sequencePtrFirst = gwiz::SequenceManager::Instance()->getSequence(sequenceFirst);
	auto sequencePtrSecond = gwiz::SequenceManager::Instance()->getSequence(sequenceSecond);
	ASSERT_NE(sequencePtrFirst.get(), sequencePtrSecond.get());
}

TEST_F(SequenceManagerTest, GetDifferentSequences)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
	auto sequencePtrFirst = gwiz::SequenceManager::Instance()->getSequence(sequenceFirst);
	auto sequencePtrSecond = gwiz::SequenceManager::Instance()->getSequence(sequenceSecond);
	ASSERT_STRNE(sequenceFirst, sequencePtrSecond->getSequence());
}

TEST_F(SequenceManagerTest, GetSequenceCountNoneAdded)
{
	ASSERT_EQ(0, gwiz::SequenceManager::Instance()->getSequenceCount());
}

TEST_F(SequenceManagerTest, GetSequenceCountOneAdded)
{
    const char* sequence = "ATGC";
	auto sequencePtr1 = gwiz::SequenceManager::Instance()->getSequence(sequence);
	auto sequencePtr2 = gwiz::SequenceManager::Instance()->getSequence(sequence);
	ASSERT_EQ(1, gwiz::SequenceManager::Instance()->getSequenceCount());
}

TEST_F(SequenceManagerTest, GetSequenceCountTwoAdded)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
	auto sequencePtr1 = gwiz::SequenceManager::Instance()->getSequence(sequenceFirst);
	auto sequencePtr2 = gwiz::SequenceManager::Instance()->getSequence(sequenceSecond);
	ASSERT_EQ(2, gwiz::SequenceManager::Instance()->getSequenceCount());
}

TEST_F(SequenceManagerTest, TestClearSequences)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
    auto sequencePtr = gwiz::SequenceManager::Instance()->getSequence(sequenceFirst);
    auto sequencePtr2 = gwiz::SequenceManager::Instance()->getSequence(sequenceSecond);
    ASSERT_EQ(2, gwiz::SequenceManager::Instance()->getSequenceCount());
    gwiz::SequenceManager::Instance()->clearSequences();
    ASSERT_EQ(0, gwiz::SequenceManager::Instance()->getSequenceCount());
}

#endif //GWIZ_TESTS_SEQUENCEMANAGERTESTS_HPP
