#ifndef GRAPHITE_VARIANT_H
#define GRAPHITE_VARIANT_H

#include "VCFWriter.h"
#include "core/util/Noncopyable.hpp"
#include "core/util/Types.h"
#include "core/sample/Sample.h"
#include "core/allele/Allele.h"
// #include "core/graph/Node.h"

#include <vector>
#include <memory>

namespace graphite
{
	class Variant : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Variant > SharedPtr;
		Variant(const std::string& variantLine, VCFWriter::SharedPtr vcfWriterPtr);
		~Variant();

		std::string getChromosome() { return m_chrom; }
		position getPosition() { return m_position; }
		Allele::SharedPtr getReferenceAllelePtr() { return this->m_reference_allele_ptr; }
		std::vector< Allele::SharedPtr > getAlternateAllelePtrs() { return this->m_alternate_allele_ptrs; }

		void writeVariant();

	private:
		void parseColumns();
		void setAlleles();
		std::string getGraphiteCounts(const std::string& sampleName);
		/* std::vector< Node::SharedPtr > getReferenceNodePtrs(); */

		VCFWriter::SharedPtr m_vcf_writer_ptr;
		std::unordered_map< std::string, std::string > m_columns;
		std::string m_variant_line;
		std::string m_chrom;
		position m_position;
		Allele::SharedPtr m_reference_allele_ptr; // make sure to figure out  a way to keep track of breaking up the alleles
		std::vector< Allele::SharedPtr > m_alternate_allele_ptrs;
		bool m_skip_adjudication;

	};
}

#endif //GRAPHITE_VARIANT_H
