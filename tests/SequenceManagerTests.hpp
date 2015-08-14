#ifndef GRAPHITE_TESTS_SEQUENCEMANAGERTESTS_HPP
#define GRAPHITE_TESTS_SEQUENCEMANAGERTESTS_HPP

#include "core/sequence/SequenceManager.h"

class SequenceManagerTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        graphite::SequenceManager::Instance()->clearSequences();
    }
};

// Tests factorial of negative numbers.
TEST_F(SequenceManagerTest, GetSequence)
{
	const char* sequence = "ATGC";
	auto sequencePtr = graphite::SequenceManager::Instance()->getSequence(sequence);
	ASSERT_STREQ(sequence, sequencePtr->getSequence());
}

TEST_F(SequenceManagerTest, CheckSameSequenceSamePtr)
{
	const char* sequence = "ATGC";
	auto sequencePtrFirst = graphite::SequenceManager::Instance()->getSequence(sequence);
	auto sequencePtrSecond = graphite::SequenceManager::Instance()->getSequence(sequence);
	ASSERT_EQ(sequencePtrFirst.get(), sequencePtrSecond.get());
}

TEST_F(SequenceManagerTest, CheckDifferentSequenceDifferentPtr)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
	auto sequencePtrFirst = graphite::SequenceManager::Instance()->getSequence(sequenceFirst);
	auto sequencePtrSecond = graphite::SequenceManager::Instance()->getSequence(sequenceSecond);
	ASSERT_NE(sequencePtrFirst.get(), sequencePtrSecond.get());
}

TEST_F(SequenceManagerTest, GetDifferentSequences)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
	auto sequencePtrFirst = graphite::SequenceManager::Instance()->getSequence(sequenceFirst);
	auto sequencePtrSecond = graphite::SequenceManager::Instance()->getSequence(sequenceSecond);
	ASSERT_STRNE(sequenceFirst, sequencePtrSecond->getSequence());
}

TEST_F(SequenceManagerTest, GetSequenceCountNoneAdded)
{
	ASSERT_EQ(0, graphite::SequenceManager::Instance()->getSequenceCount());
}

TEST_F(SequenceManagerTest, GetSequenceCountOneAdded)
{
    const char* sequence = "ATGC";
	auto sequencePtr1 = graphite::SequenceManager::Instance()->getSequence(sequence);
	auto sequencePtr2 = graphite::SequenceManager::Instance()->getSequence(sequence);
	ASSERT_EQ(1, graphite::SequenceManager::Instance()->getSequenceCount());
}

TEST_F(SequenceManagerTest, GetSequenceCountTwoAdded)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
	auto sequencePtr1 = graphite::SequenceManager::Instance()->getSequence(sequenceFirst);
	auto sequencePtr2 = graphite::SequenceManager::Instance()->getSequence(sequenceSecond);
	ASSERT_EQ(2, graphite::SequenceManager::Instance()->getSequenceCount());
}

TEST_F(SequenceManagerTest, TestClearSequences)
{
    const char* sequenceFirst = "ATGC";
    const char* sequenceSecond = "ATGCATG";
    auto sequencePtr = graphite::SequenceManager::Instance()->getSequence(sequenceFirst);
    auto sequencePtr2 = graphite::SequenceManager::Instance()->getSequence(sequenceSecond);
    ASSERT_EQ(2, graphite::SequenceManager::Instance()->getSequenceCount());
    graphite::SequenceManager::Instance()->clearSequences();
    ASSERT_EQ(0, graphite::SequenceManager::Instance()->getSequenceCount());
}

#endif //GRAPHITE_TESTS_SEQUENCEMANAGERTESTS_HPP
