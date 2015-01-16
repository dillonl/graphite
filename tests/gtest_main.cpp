#include <stdio.h>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/directed_graph.hpp"

#include "VCFFileTests.hpp"

#include "BuildGraphTests.hpp"
#include "VariantsTest.hpp"
#include "ReferenceTest.hpp"
#include "RegionTests.hpp"
#include "FileTests.hpp"

#include "gtest/gtest.h"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
