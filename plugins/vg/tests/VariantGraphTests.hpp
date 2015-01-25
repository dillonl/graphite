#ifndef GWIZ_VG_VARIANT_GRAPH_TESTS_HPP
#define GWIZ_VG_VARIANT_GRAPH_TESTS_HPP

#include "gtest/gtest.h"

#include "TestConfig.h"
#include "core/variants/IVariant.h"
#include "core/reference/FastaReference.h"

namespace
{


	class VariantGraphTest : public ::testing::Test
	{
	protected:
		// You can remove any or all of the following functions if its body
		// is empty.

		VariantGraphTest()
		{
				// You can do set-up work for each test here.
		}

		virtual ~VariantGraphTest()
		{
			// You can do clean-up work that doesn't throw exceptions here.
	    }

		// If the constructor and destructor are not enough for setting up
		// and cleaning up each test, you can define the following methods:

		virtual void SetUp()
		{
			// Code here will be called immediately after the constructor (right
			// before each test).
		}

		virtual void TearDown()
		{
			// Code here will be called immediately after each test (right
			// before the destructor).
		}

		// Objects declared here can be used by all tests in the test case for Foo.
	};

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(VariantGraphTest, ConstructGraph)
	{
		/*
		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20:100000-1000000");
		std::string fastaPath = TEST_FASTA_FILE;
		gwiz::FastaReference::SharedPtr fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		const char* reference = fastaReferencePtr->getSequence();
		*/
	}
}

#endif //VARIANT_GRAPH_TESTS_HPP
