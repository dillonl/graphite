#ifndef GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
#define GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
#include <chrono>
#include <thread>

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "plugins/adjudicator/graph/AdjudicatorGraphTests.hpp"

#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/reference/FastaReference.h"

/*
namespace
{

	class VGTest : public gwiz::vg::AdjudicatorGraph
	{
	public:
		typedef std::shared_ptr< VGTest > VGTestPtr;

		VGTest(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr) :
			VariantGraph(referencePtr, variantListPtr)
			{

			}

		~VGTest()
			{
			}

		void constructGraph() override
			{
				gwiz::vg::VariantGraph::constructGraph();
			}

		bool getNextCompoundVariant(gwiz::Variant::SharedPtr& variant) override
			{
				return gwiz::vg::VariantGraph::getNextCompoundVariant(variant);
			}

		gwiz::Variant::SharedPtr buildCompoundVariant(const gwiz::position startPosition, const std::string& referenceString, const std::vector< gwiz::Variant::SharedPtr >& variants) override
			{
				return gwiz::vg::VariantGraph::buildCompoundVariant(startPosition, referenceString, variants);
			}
	};

}
*/
#endif //GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
