#ifndef GWIZ_VARIANT_H
#define GWIZ_VARIANT_H

#include <exception>
#include <cstring>
#include <string>
#include <vector>
#include <memory>

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
			auto variant = std::make_shared< Variant >();
			variant->m_variant_type = VARIANT_TYPE::SNP;
			if (!boost::spirit::qi::parse(vcf_line, end_line, parser, variant->m_chrom, variant->m_position, variant->m_id, variant->m_ref, variant->m_alt))
			{
				throw "An invalid line caused an exception. Please correct the input and try again";
			}
			return variant;
		}

		VARIANT_TYPE getVariantType() const { return m_variant_type; }
		std::string getChrom() const { return m_chrom; }
		uint32_t getPosition() const { return m_position; }
		std::string getID() const { return m_id; }
		std::vector< std::string > const getRef() { return m_ref; }
		std::vector< std::string > const getAlt() { return m_alt; }

		size_t getSmallestAlleleSize() override; // returns the smallest allele in this variant (including reference allele)
		size_t getLargestAlleleSize() override; // returns the largest allele in this variant (including reference allele)
	private:

		VARIANT_TYPE m_variant_type;
		uint32_t m_position;
		std::string m_chrom;
		std::string m_id;
		std::vector< std::string > m_ref;
		std::vector< std::string > m_alt;
	};

}

#endif //GWIZ_VARIANT_H
