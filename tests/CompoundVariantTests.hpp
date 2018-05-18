#ifndef GRAPHITE_COMPOUNDVARIANTTESTS_HPP
#define GRAPHITE_COMPOUNDVARIANTTESTS_HPP

#include <vector>

#include "TestConfig.h"
#include "TestClasses.hpp"

#include "core/allele/Allele.h"
#include "core/allele/EquivalentAllele.h"
#include "core/variant/Variant.h"
#include "core/variant/VCFFileReader.h"

namespace
{
	using namespace graphite;

	void compareAlleles(IAllele::SharedPtr allele1, IAllele::SharedPtr allele2)
	{
		// std::cout << "----" << std::endl;
		auto eAllele1 = std::dynamic_pointer_cast< EquivalentAllele >(allele1);
		auto eAllele2 = std::dynamic_pointer_cast< EquivalentAllele >(allele2);
		if (eAllele1 != nullptr && eAllele2 != nullptr)
		{
			// std::cout << "ea1(r): " << eAllele1->getSequence() << " ea2(r): " << eAllele2->getSequence() << std::endl;
			ASSERT_STREQ(eAllele1->getSequence(), eAllele2->getSequence());
			for (uint32_t i = 0; i < eAllele1->getAllAlleles().size(); ++i)
			{
				auto eqAllele1 = eAllele1->getAllAlleles()[i];
				auto eqAllele2 = eAllele2->getAllAlleles()[i];
				auto eq1MetaDataPtr = eqAllele1->getAlleleMetaData();
				auto eq2MetaDataPtr = eqAllele2->getAlleleMetaData();
				ASSERT_STREQ(eqAllele1->getSequence(), eqAllele2->getSequence());
				ASSERT_EQ(eq1MetaDataPtr->getPaddingPrefix(), eq2MetaDataPtr->getPaddingPrefix());
				ASSERT_EQ(eq1MetaDataPtr->getPaddingSuffix(), eq2MetaDataPtr->getPaddingSuffix());
			}
			ASSERT_EQ(eAllele1->getAllAlleles().size(), eAllele2->getAllAlleles().size());
		}
		else
		{
			auto allele1MetaDataPtr = allele1->getAlleleMetaData();
			auto allele2MetaDataPtr = allele2->getAlleleMetaData();
			ASSERT_STREQ(allele1->getSequence(), allele2->getSequence());
			ASSERT_EQ(allele1MetaDataPtr->getPaddingPrefix(), allele2MetaDataPtr->getPaddingPrefix());
			ASSERT_EQ(allele1MetaDataPtr->getPaddingSuffix(), allele2MetaDataPtr->getPaddingSuffix());
		}
		// std::cout << "----" << std::endl;
	}

	void testCompoundVariants(std::vector< std::string > variantLines, IReference::SharedPtr referencePtr, std::vector< IAllele::SharedPtr > refAllelePtr, std::vector< std::vector< IAllele::SharedPtr > > altAlleles)
	{
		std::vector< IVariant::SharedPtr > variantPtrs;
		for (auto& variantLine : variantLines)
		{
			variantPtrs.emplace_back(Variant::BuildVariant(variantLine, referencePtr, 100));
		}
		auto variantsListPtr = std::make_shared< VariantList >(variantPtrs, nullptr);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		uint32_t variantCounter = 0;
		IVariant::SharedPtr variantPtr;
		while (variantsListPtr->getNextVariant(variantPtr))
		{
			compareAlleles(variantPtr->getRefAllelePtr(), refAllelePtr[variantCounter]);
			for (size_t i = 0; i < altAlleles[variantCounter].size(); ++i)
			for (size_t i = 0; i < variantPtr->getAltAllelePtrs().size(); ++i)
			{
				// std::cout << variantPtr->getPosition() << std::endl;
				compareAlleles(variantPtr->getAltAllelePtrs()[i], altAlleles[variantCounter][i]);
			}
			ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), altAlleles[variantCounter].size());
			++variantCounter;
		}
	}

	TEST(CompoundVariantTests, CompoundVariantSimple)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", 1,156, graphite::Region::BASED::ONE);
		regionPtr->setBased(graphite::Region::BASED::ONE);
		std::string sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT";
		auto referencePtr = std::make_shared< graphite::Reference >(sequence, regionPtr);

		std::vector< std::string > variantLines = { "1\t60\t.\tTAAGT\tT\t.\t.\t.", "1\t61\t.\tAA\tAGAAG\t.\t.\t.\t.", "1\t63\t.\tGTTAAC\tGTAG,G\t.\t.\t." };
		std::vector< IAllele::SharedPtr > refAllelePtrList = { std::make_shared< Allele >("TAAGTTAAC", std::make_shared< AlleleMetaData >(0,0)) };
		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		std::vector< IAllele::SharedPtr > altAllelePtrs;
		altAllelePtrs.emplace_back(std::make_shared< Allele >("TTAAC", std::make_shared< AlleleMetaData >(0,4)));
		altAllelePtrs.emplace_back(std::make_shared< Allele >("TAGAAGGTTAAC", std::make_shared< AlleleMetaData >(1,6)));
		altAllelePtrs.emplace_back(std::make_shared< Allele >("TAAGTAG", std::make_shared< AlleleMetaData >(3,0)));
		altAllelePtrs.emplace_back(std::make_shared< Allele >("TAAG", std::make_shared< AlleleMetaData >(3,0)));
		altAllelePtrsList.emplace_back(altAllelePtrs);
		testCompoundVariants(variantLines, referencePtr, refAllelePtrList, altAllelePtrsList);
	}

	TEST(CompoundVariantTests, BuildCompoundVariant)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", 1,156, graphite::Region::BASED::ONE);
		regionPtr->setBased(graphite::Region::BASED::ONE);
		std::string sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT";
		auto referencePtr = std::make_shared< graphite::Reference >(sequence, regionPtr);

		std::vector< std::string > variantLines = { "1\t46\t.\tTAATATATGTAATATATATTATATATGTAATATAATATATGTAAT\tT\t.\t.\t.", "1\t55\t.\tTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT\tT\t.\t.\t.", "1\t104\t.\tGT\tG\t.\t.\t." };
		std::vector< IAllele::SharedPtr > refAllelePtrs = {
			std::make_shared< Allele >("TAATATATGTAATATATATTATATATGTAATATAATATATGTAAT", std::make_shared< AlleleMetaData >(0, 65)),
			std::make_shared< Allele >("TAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", std::make_shared< AlleleMetaData >(9, 0)),
			std::make_shared< Allele >("GT", std::make_shared< AlleleMetaData >(58, 50))
		};

		auto refPtr = std::make_shared< EquivalentAllele >("TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", refAllelePtrs);
		std::vector< IAllele::SharedPtr > refAllelePtrList = { refPtr };

		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		std::vector< IAllele::SharedPtr > altAllelePtrs = {
			std::make_shared< Allele >("TATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", std::make_shared< AlleleMetaData >(0, 65)),
			std::make_shared< Allele >("TAATATATGT", std::make_shared< AlleleMetaData >(9, 0)),
			std::make_shared< Allele >("TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", std::make_shared< AlleleMetaData >(58, 50))
			};
		altAllelePtrsList.emplace_back(altAllelePtrs);
		testCompoundVariants(variantLines, referencePtr, refAllelePtrList, altAllelePtrsList);
	}

	TEST(CompoundVariantTests, BuildSemanticCompoundVariant)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", 1,156, graphite::Region::BASED::ONE);
		regionPtr->setBased(graphite::Region::BASED::ONE);
		std::string sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT";
		auto referencePtr = std::make_shared< graphite::Reference >(sequence, regionPtr);

		std::vector< std::string > variantLines = {
			"1\t105774\t.\tGT\tG\t.\t.\t.",
			"1\t105774\t.\tGTT\tGT\t.\t.\t." };

		std::vector< IAllele::SharedPtr > refAllelePtrs = {
			std::make_shared< Allele >("GT", std::make_shared< AlleleMetaData >(0, 1)),
			std::make_shared< Allele >("GTT", std::make_shared< AlleleMetaData >(0, 0)),
		};
		auto refPtr = std::make_shared< EquivalentAllele >("GTT", refAllelePtrs);
		std::vector< IAllele::SharedPtr > refAllelePtrsList = { refPtr };

		std::vector< IAllele::SharedPtr > altAlleleEqPtrs = {
			std::make_shared< Allele >("GT", std::make_shared< AlleleMetaData >(0, 1)),
			std::make_shared< Allele >("GT", std::make_shared< AlleleMetaData >(0, 0)),
		};
		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		auto altPtr = std::make_shared< EquivalentAllele >("GT", altAlleleEqPtrs);

		std::vector< IAllele::SharedPtr > altAllelePtrs = { altPtr };
		altAllelePtrsList.emplace_back(altAllelePtrs);

		testCompoundVariants(variantLines, referencePtr, refAllelePtrsList, altAllelePtrsList);
		// testCompoundVariants(variantLines, referencePtr, refAllelePtrList, altAllelePtrsList);
	}

	TEST(CompoundVariantTests, BuildCompoundVariant3)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", 1,156, graphite::Region::BASED::ONE);
		regionPtr->setBased(graphite::Region::BASED::ONE);
		std::string sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT";
		auto referencePtr = std::make_shared< graphite::Reference >(sequence, regionPtr);

		std::vector< std::string > variantLines = {
			"20\t69506\t.\tG\tGACAC\t.\t.\t.",
			"20\t69506\t.\tGACACACACACACACACACACACACACACA\tGACACACACACACACACACACACACACACACACA\t.\t.\t."
		};
		// ref allele
		std::vector< IAllele::SharedPtr > refAllelePtrs = {
			std::make_shared< Allele >("G", std::make_shared< AlleleMetaData >(0, 29)),
			std::make_shared< Allele >("GACACACACACACACACACACACACACACA", std::make_shared< AlleleMetaData >(0, 0)),
		};
		auto refPtr = std::make_shared< EquivalentAllele >("GACACACACACACACACACACACACACACA", refAllelePtrs);
		std::vector< IAllele::SharedPtr > refAllelePtrsList = { refPtr };

		// alt allele
		std::vector< IAllele::SharedPtr > altAlleleEqPtrs = {
			std::make_shared< Allele >("GACACACACACACACACACACACACACACACACA", std::make_shared< AlleleMetaData >(0, 29)),
			std::make_shared< Allele >("GACACACACACACACACACACACACACACACACA", std::make_shared< AlleleMetaData >(0, 0)),
		};
		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		auto altPtr = std::make_shared< EquivalentAllele >("GACACACACACACACACACACACACACACACACA", altAlleleEqPtrs);

		std::vector< IAllele::SharedPtr > altAllelePtrs = { altPtr };
		altAllelePtrsList.emplace_back(altAllelePtrs);

		testCompoundVariants(variantLines, referencePtr, refAllelePtrsList, altAllelePtrsList);
	}
	TEST(CompoundVariantTests, BuildGraph4)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", 1,156, graphite::Region::BASED::ONE);
		regionPtr->setBased(graphite::Region::BASED::ONE);
		std::string sequence = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT";
		auto referencePtr = std::make_shared< graphite::Reference >(sequence, regionPtr);

		std::vector< std::string > variantLines = {
			"20\t72104\t.\tTA\tT\t.\t.\t.",
			"20\t72104\t.\tTAA\tTA\t.\t.\t.",
			"20\t72719\t.\tC\tT\t.\t.\t." };

		// ref allele
		std::vector< IAllele::SharedPtr > refAllelePtrs = {
			std::make_shared< Allele >("TA", std::make_shared< AlleleMetaData >(0, 1)),
			std::make_shared< Allele >("TAA", std::make_shared< AlleleMetaData >(0, 0)),
		};
		auto refPtr1 = std::make_shared< EquivalentAllele >("TAA", refAllelePtrs);
		auto refPtr2 = std::make_shared< Allele >("C", std::make_shared< AlleleMetaData >(0, 0));
		std::vector< IAllele::SharedPtr > refAllelePtrsList = { refPtr1, refPtr2 };


		// alt allele
		std::vector< IAllele::SharedPtr > altAlleleEqPtrs = {
			std::make_shared< Allele >("TA", std::make_shared< AlleleMetaData >(0, 1)),
			std::make_shared< Allele >("TA", std::make_shared< AlleleMetaData >(0, 0)),
		};
		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		auto altPtr1 = std::make_shared< EquivalentAllele >("TA", altAlleleEqPtrs);
		auto altPtr2 = std::make_shared< Allele >("T", std::make_shared< AlleleMetaData >(0, 0));

		std::vector< IAllele::SharedPtr > altAllelePtrs1 = { altPtr1 };
		std::vector< IAllele::SharedPtr > altAllelePtrs2 = { altPtr2 };

		altAllelePtrsList.emplace_back(altAllelePtrs1);
		altAllelePtrsList.emplace_back(altAllelePtrs2);

		testCompoundVariants(variantLines, referencePtr, refAllelePtrsList, altAllelePtrsList);
	}

}
#endif //GRAPHITE_COMPOUNDVARIANTTESTS_HPP
