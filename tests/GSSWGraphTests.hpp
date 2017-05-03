#ifndef GRAPHITE_GSSWGRAPHTESTS_HPP
#define GRAPHITE_GSSWGRAPHTESTS_HPP

#include "TestConfig.h"

#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/graph/GSSWGraph.h"

/*
TEST(GSSWGraphTests, GetBaseOneRegion)
{
	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);
	auto regionPtr = std::make_shared< graphite::Region >("1:1-10", graphite::Region::BASED::ONE);
	std::string regionString = referencePtr->getSequenceFromRegion(regionPtr);
	ASSERT_STREQ(regionString.c_str(), "CTATGATGTT");
}
*/

TEST(GSSWGraphTests, GSSWSimpleGraph)
{
    uint32_t readLength = 6;
	std::string vcfLine = "1\t10\trs11575897\tT\tG\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line
	auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);
	auto variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), referencePtr, readLength);

	std::vector< graphite::IVariant::SharedPtr > variantPtrs = {variantPtr};
	auto variantListPtr = std::make_shared< graphite::VariantList >(variantPtrs, referencePtr);

    auto gsswRegionPtr = std::make_shared< graphite::Region >("1", variantPtr->getPosition() - readLength, variantPtr->getPosition() + readLength, graphite::Region::BASED::ONE);
	auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, gsswRegionPtr, 1, 1, 1, 1, 1);
	gsswGraphPtr->constructGraph();

	gssw_graph* gsswPtr = gsswGraphPtr->getGSSWGraph();
	ASSERT_TRUE(gsswPtr->size == 4);
	// TGATGT T GATGGAA
    ASSERT_STREQ(gsswPtr->nodes[0]->seq, "TGATGT");
	ASSERT_STREQ(gsswPtr->nodes[1]->seq, "G");
	ASSERT_STREQ(gsswPtr->nodes[2]->seq, "T");
    ASSERT_STREQ(gsswPtr->nodes[3]->seq, "GATGGA"); // GATGGA
}

TEST(GSSWGraphTests, GSSWSimpleLargeVariant)
{
	uint32_t readLength = 3;
	std::string vcfLine = "1\t20\trs11575897\tGACCAAACGTCGTTAGGCCAGTTTTCTGGT\tG\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line
	auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);
	auto variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), referencePtr, readLength);

	std::vector< graphite::IVariant::SharedPtr > variantPtrs = {variantPtr};
	auto variantListPtr = std::make_shared< graphite::VariantList >(variantPtrs, referencePtr);

	auto gsswRegionPtr = std::make_shared< graphite::Region >("1", variantPtr->getPosition() - readLength, variantPtr->getPosition() + 29 + readLength, graphite::Region::BASED::ONE);
	auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, gsswRegionPtr, 1, 1, 1, 1, 1);
	gsswGraphPtr->constructGraph();

	gssw_graph* gsswPtr = gsswGraphPtr->getGSSWGraph();
	ASSERT_TRUE(gsswPtr->size == 4);
	ASSERT_STREQ(gsswPtr->nodes[0]->seq, "ACT");
	ASSERT_STREQ(gsswPtr->nodes[1]->seq, "G");
	ASSERT_STREQ(gsswPtr->nodes[2]->seq, "GACCNNNGGT"); // GACCAAACGTCGTTAGGCCAG
	ASSERT_STREQ(gsswPtr->nodes[3]->seq, "CGT");
}

#endif
