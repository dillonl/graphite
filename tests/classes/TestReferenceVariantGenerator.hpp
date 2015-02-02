#ifndef GWIZ_TEST_TESTREFERENCEVARIANTGENERATOR_HPP
#define GWIZ_TEST_TESTREFERENCEVARIANTGENERATOR_HPP

#include <vector>
#include <boost/algorithm/string/join.hpp>

namespace gwiz
{
	namespace testing
	{
		class TestReferenceVariantGenerator
		{
		public:
		TestReferenceVariantGenerator(const std::string& reference, const std::string& chrom, position startPosition) :
			    m_reference(reference)
			{
				m_region = std::make_shared< Region >(chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(startPosition + reference.size()));
			}
			~TestReferenceVariantGenerator()
			{
			}

			void addVariant(position position, const std::string& id, uint32_t refLength, const std::vector< std::string >& alts)
			{
				if (position > this->m_region->getEndPosition())
				{
					throw "Position of variant is past the end position of the reference region";
				}
				std::string ref = std::string(this->m_reference[position - this->m_region->getStartPosition()], 1);
				if (std::find(alts.begin(), alts.end(), ref) != alts.end())
				{
					throw "Reference and alt are the same at position: " + std::to_string(position);
				}
				std::string variantString = boost::algorithm::join(alts, ",");
				std::string variantLine = m_region->getReferenceID() + "\t" + std::to_string(position) + "\t" + id + "\t" + std::string(m_reference.c_str() + (position - m_region->getStartPosition()), refLength) + "\t" + variantString + "\t";

				m_variant_list.push_back(Variant::BuildVariant(variantLine.c_str(), m_vcf_parser));
				std::cout << "variant added: " << m_variant_list.size() << std::endl;
			}

			TestVariantList::SharedPtr getVariants()
			{
				return std::make_shared< TestVariantList >(m_variant_list);
			}

			TestReference::SharedPtr getReference()
			{
				return std::make_shared< TestReference >(m_reference, this->m_region->getRegionString());
			}

		private:
			std::string m_reference;
			Region::SharedPtr m_region;

			std::vector< Variant::SharedPtr > m_variant_list;
			VariantParser< const char* > m_vcf_parser;


		};
	}

}

#endif //GWIZ_TEST_TESTREFERENCEVARIANTGENERATOR_HPP
