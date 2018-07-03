#ifndef GRAPHITE_GRAPHPROCESSOR_H
#define GRAPHITE_GRAPHPROCESSOR_H

#include "core2/util/Noncopyable.hpp"
#include "core2/util/ThreadPool.hpp"
#include "core2/region/Region.h"
#include "core2/reference/FastaReference.h"
#include "core2/vcf/VCFReader.h"
#include "core2/bam/BamReader.h"
#include "Graph.h"

#include "api/BamAlignment.h"

#include <memory>
#include <mutex>
#include <condition_variable>
#include <future>
#include <unordered_set>
#include <thread>
#include <queue>
#include <algorithm>
#include <functional>
#include <atomic>

namespace graphite
{
	class GraphProcessor : private Noncopyable
	{
	public:
		typedef std::shared_ptr< GraphProcessor > SharedPtr;
		GraphProcessor(FastaReference::SharedPtr fastaReferencePtr, const std::vector< BamReader::SharedPtr >& bamReaderPtrs, const std::vector< VCFReader::SharedPtr >& vcfReaderPtrs, uint32_t matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue);
		~GraphProcessor();

		void processVariants();

	private:
		void adjudicateVariants(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t graphSpacing);
		void adjudicateVariants2(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t graphSpacing);
        void getAlignmentsInRegion(std::vector< std::shared_ptr< BamAlignment > >& bamAlignmentPtrs, std::vector< Region::SharedPtr > regionPtrs, bool getFlankingUnalignedReads);
		FastaReference::SharedPtr m_fasta_reference_ptr;
		std::vector< BamReader::SharedPtr > m_bam_reader_ptrs;
		std::vector< VCFReader::SharedPtr > m_vcf_reader_ptrs;
		std::unordered_map< std::string, Sample::SharedPtr > m_bam_sample_ptrs;
		uint32_t m_flanking_padding;
		uint32_t m_match_value;
		uint32_t m_mismatch_value;
		uint32_t m_gap_open_value;
		uint32_t m_gap_extension_value;

		ThreadPool m_threadpool;
	};
}

#endif //GRAPHITE_GRAPHPROCESSOR_H
