#ifndef GWIZ_VARIANT_H
#define GWIZ_VARIANT_H

#include <exception>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <map>

#include "IVariant.h"
#include "VCFParser.hpp"
#include "core/reference/Reference.h"

namespace gwiz
{

	enum class VARIANT_TYPE
	{
		SNP,
		INS,
		DEL,
		DUP,
		INV,
	};

	class Variant : public IVariant
	{
	public:
		typedef std::shared_ptr<Variant> SharedPtr;
		Variant();
		~Variant();

		inline static Variant::SharedPtr BuildVariant(const char* vcf_line, VariantParser< const char* >& parser)
		{
			const char* end_line = static_cast< const char* >(memchr(vcf_line, '\n', std::numeric_limits< position >::max()));
			auto variantPtr = std::make_shared< Variant >();
			variantPtr->m_variant_type = VARIANT_TYPE::SNP;
			if (!boost::spirit::qi::parse(vcf_line, end_line, parser, variantPtr->m_chrom, variantPtr->m_position, variantPtr->m_id, variantPtr->m_ref, variantPtr->m_alt))
			{
				throw "An invalid line in the VCF caused an exception. Please correct the input and try again";
			}
			for (auto alternateString : variantPtr->getAlt())
			{
				variantPtr->setVCFLineFromAlternate(alternateString, vcf_line);
			}
			return variantPtr;
		}

		inline void setVCFLineFromAlternate(const std::string& alt, const char* vcf_line)
		{
			this->m_vcf_lines_map[alt] = vcf_line;
		}

		inline const char* getVCFLineFromAlternate(const std::string& alt)
		{
			if (this->m_vcf_lines_map.find(alt) == this->m_vcf_lines_map.end())
			{
				return "";
			}
			return this->m_vcf_lines_map[alt];
		}

		inline void increaseCount(const std::string& allele)
		{
			size_t count = 0;
			if (this->m_allele_count.find(allele) == this->m_allele_count.end())
			{
				count = this->m_allele_count[allele];
			}
			this->m_allele_count[allele] = count + 1;
			++this->m_total_allele_count;
		}

		uint32_t getAlleleCount(const std::string& allele)
		{
			size_t count = 0;
			if (this->m_allele_count.find(allele) != this->m_allele_count.end())
			{
				count = this->m_allele_count[allele];
			}
			return count;
		}

		void generateAlleleCoveragePercentages(std::vector< std::tuple< std::string, uint32_t > >& allelePercentages)
		{
			for (auto alleleCountIter : this->m_allele_count)
			{
				auto allelePercent = (this->m_total_allele_count > 0) ? alleleCountIter.second / this->m_total_allele_count : 0;
				std::tuple< std::string, uint32_t > alleleTuple(alleleCountIter.first, allelePercent);
				allelePercentages.emplace_back(alleleTuple);
			}
			std::sort(allelePercentages.begin(), allelePercentages.end(), [](const std::tuple< std::string, uint32_t >& lhs, const std::tuple< std::string, uint32_t >& rhs)
					  {
						  return std::get< 1 >(lhs) < std::get< 1 >(rhs);
					  });
		}

		void printAlleleCount()
		{
			for (auto alleleCounter : this->m_allele_count)
			{
				std::cout << "allele: " << alleleCounter.first << " <" << alleleCounter.second << ">" << std::endl;
			}
			std::cout << "allele: " << getRef() << " <" << this->m_allele_count[getRef()] << ">" << std::endl;
		}

		VARIANT_TYPE getVariantType() const { return m_variant_type; }
		std::string getChrom() const { return m_chrom; }
		uint32_t getPosition() const { return m_position; }
		std::string getID() const { return m_id; }
		std::string const getRef() { return m_ref[0]; }
		std::vector< std::string > const getAlt() { return m_alt; }
		std::map< std::string, uint32_t > m_allele_count;
		uint32_t m_total_allele_count; // an efficiency that technically could be calculated from m_allele_count

		size_t getSmallestAlleleSize() override; // returns the smallest allele in this variant (including reference allele)
		size_t getLargestAlleleSize() override; // returns the largest allele in this variant (including reference allele)

		std::string toString() override;
	private:

		VARIANT_TYPE m_variant_type;
		uint32_t m_position;
		std::string m_chrom;
		std::string m_id;
		std::vector< std::string > m_ref;
		std::vector< std::string > m_alt;
		std::map< std::string, const char* > m_vcf_lines_map;
	};

}

#endif //GWIZ_VARIANT_H
