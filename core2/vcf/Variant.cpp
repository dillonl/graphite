#include "Variant.h"
#include "core2/util/Utility.h"
#include "core2/util/Types.h"

namespace graphite
{

	Variant::Variant(const std::string& variantLine, VCFWriter::SharedPtr vcfWriterPtr) :
		m_variant_line(variantLine),
		m_vcf_writer_ptr(vcfWriterPtr),
		m_skip_adjudication(false)
	{
		parseColumns();
		setAlleles();
	}

	Variant::~Variant()
	{
	}

	void Variant::writeVariant()
	{
	}

	void Variant::parseColumns()
	{
		auto samplePtrs = this->m_vcf_writer_ptr->getSamplePtrs();
		std::vector< std::string > columns;
		split(this->m_variant_line, '\t', columns);
		if (columns.size() != (STANDARD_VCF_COLUMN_NAMES.size() + samplePtrs.size()))
		{
			std::cout << "Invalid VCF at line: " << this->m_variant_line << std::endl;
			exit(EXIT_FAILURE);
		}
		for (uint32_t i = 0; i < STANDARD_VCF_COLUMN_NAMES.size(); ++i)
		{
			auto standardColumn = STANDARD_VCF_COLUMN_NAMES[i];
			m_columns[standardColumn] = columns[i];
		}
		for (uint32_t i = 0; i < samplePtrs.size(); ++i)
		{
			m_columns[samplePtrs[i]->getName()] = columns[i];
		}
	}

	void Variant::setAlleles()
	{
		this->m_reference_allele_ptr = std::make_shared< Allele >(m_columns["REF"]);
		std::vector< std::string > alts;
		split(m_columns["ALT"], '\t', alts);
		this->m_alternate_allele_ptrs.clear();
		for (auto alt : alts)
		{
			auto altAllelePtr = std::make_shared< Allele >(alt);
			this->m_alternate_allele_ptrs.emplace_back(altAllelePtr);
		}
	}
}
