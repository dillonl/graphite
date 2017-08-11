#ifndef GRAPHITE_BAMHEADERREADERTESTS_HPP
#define GRAPHITE_BAMHEADERREADERTESTS_HPP

#include "core/file/BamHeaderReader.h"
#include "TestConfig.h"

#include <iostream>

/*
 * Check that RNAME's are included in the header.
 */ 
class BamHeaderReaderTest : public ::testing::Test
{
protected:
    std::string bamPath = {TEST_BAM_FILE};
    graphite::BamHeaderReader bamReader = graphite::BamHeaderReader(bamPath);

    virtual void SetUp ()
    {
        bamReader.open();
    }

    virtual void TearDown ()
    {
        bamReader.close();
    }
};

TEST_F (BamHeaderReaderTest, addReadGroups_toSamHeader)
{
    bamReader.addReadGroupsToSamHeader();
    ASSERT_TRUE(bamReader.containsReadGroup("ALT"));
    ASSERT_TRUE(bamReader.containsReadGroup("REF"));
}

/*
TEST_F (BamHeaderReaderTest, addRnameTo_header)
{
}
*/

#endif // GRAPHITE_BAMHEADERREADERTESTS_HPP
