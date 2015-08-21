#ifndef GRAPHITE_ADJUDICATIONTESTS_HPP
#define GRAPHITE_ADJUDICATIONTESTS_HPP

#include "TestConfig.h"

#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/variant/IVariant.h"
#include "core/adjudicator/IAdjudicator.h"

#include "plugins/gssw/graph/GSSWGraph.h"

#include <vector>

namespace
{
namespace adj_test
{
	using namespace graphite;
	using namespace graphite::gssw;
	class ReferenceTest : public IReference
	{
	public:
        ReferenceTest(const char* sequence) : m_sequence(sequence) {}
		~ReferenceTest() {}

		const char* getSequence() override
		{
			return m_sequence;
		}

		size_t getSequenceSize() override
		{
			return strlen(m_sequence);
		}

	private:
		const char* m_sequence = "";
	};

	class VariantListTest : public IVariantList
	{
	public:
		VariantListTest(std::vector< IVariant::SharedPtr > variants) : m_variant_list(variants), m_index(0) {}
		~VariantListTest() {}

		bool getNextVariant(IVariant::SharedPtr& variantPtr) override
		{
			if (this->m_variant_list.size() <= this->m_index) { variantPtr = nullptr; return false; }
			else { variantPtr = this->m_variant_list[this->m_index++]; return true; }
		}

		size_t getCount() { return this->m_variant_list.size(); }
		void sort() {}
		void printToVCF(std::ostream& out) {}

	private:
		std::vector< IVariant::SharedPtr > m_variant_list;
		uint32_t m_index;
	};

	GSSWGraph::SharedPtr getGSSWGraph(IVariantList::SharedPtr variantListPtr, IAdjudicator::SharedPtr adjudicatorPtr)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto referencePtr = std::make_shared< ReferenceTest >(seq.c_str());
		uint32_t startPosition = 1;
		uint32_t graphSize = 250;
		auto gsswGraphPtr = std::make_shared< GSSWGraph >(referencePtr, variantListPtr, startPosition, graphSize, adjudicatorPtr->getMatchValue(), adjudicatorPtr->getMisMatchValue(), adjudicatorPtr->getGapOpenValue(), adjudicatorPtr->getGapExtensionValue());
		gsswGraphPtr->constructGraph();
		return gsswGraphPtr;
	}

	TEST(AdjudicationTest, AdjudicateDualSNPMatch)
	{

	}

}
}

#endif //GRAPHITE_ADJUDICATIONTESTS_HPP
