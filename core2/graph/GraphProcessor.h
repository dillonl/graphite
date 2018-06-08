#ifndef GRAPHITE_GRAPHPROCESSOR_H
#define GRAPHITE_GRAPHPROCESSOR_H

#include "core2/util/Noncopyable.hpp"
#include "core2/region/Region.h"
#include "core2/reference/FastaReference.h"
#include "core2/vcf/VCFReader.h"
#include "core2/bam/BamReader.h"
#include "Graph.h"

#include "api/BamAlignment.h"

#include <memory>

namespace graphite
{
	class GraphProcessor : private Noncopyable
	{
	public:
		typedef std::shared_ptr< GraphProcessor > SharedPtr;
		GraphProcessor(FastaReference::SharedPtr fastaReferencePtr, const std::vector< BamReader::SharedPtr >& bamReaderPtrs, const std::vector< VCFReader::SharedPtr >& vcfReaderPtrs);
		~GraphProcessor();

		void processVariants();

	private:
		void adjudicateVariants(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t graphSpacing);
        void getAlignmentsInRegion(std::vector< std::shared_ptr< BamAlignment > >& bamAlignmentPtrs, std::vector< Region::SharedPtr > regionPtrs, bool getFlankingUnalignedReads);
		FastaReference::SharedPtr m_fasta_reference_ptr;
		std::vector< BamReader::SharedPtr > m_bam_reader_ptrs;
		std::vector< VCFReader::SharedPtr > m_vcf_reader_ptrs;
		uint32_t m_flanking_padding;
	};
}

#endif //GRAPHITE_GRAPHPROCESSOR_H
