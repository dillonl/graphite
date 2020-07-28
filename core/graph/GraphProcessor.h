#ifndef GRAPHITE_GRAPHPROCESSOR_H
#define GRAPHITE_GRAPHPROCESSOR_H

#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/vcf/VCFReader.h"
#include "core/alignment/AlignmentReader.h"
#include "core/alignment/Alignment.h"
#include "core/util/ThreadPool.hpp"
#include "core/util/GraphPrinter.h"
#include "Graph.h"

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
		GraphProcessor(FastaReference::SharedPtr fastaReferencePtr, const std::vector< AlignmentReader::SharedPtr >& alignmentReaderPtrs, const std::vector< VCFReader::SharedPtr >& vcfReaderPtrs, uint32_t matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue, bool printGraph, int32_t mappingQuality, int32_t readSampleLimit, uint32_t numberOfThreads);
		~GraphProcessor();

		void processVariants();

	private:
		void adjudicateVariants(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t graphSpacing);
        void getAlignmentsInRegion(std::vector< Alignment::SharedPtr >& alignmentPtrs, std::vector< Region::SharedPtr > regionPtrs, bool getFlankingUnalignedReads);
		FastaReference::SharedPtr m_fasta_reference_ptr;
		std::vector< AlignmentReader::SharedPtr > m_alignment_reader_ptrs;
		std::vector< VCFReader::SharedPtr > m_vcf_reader_ptrs;
		std::unordered_map< std::string, Sample::SharedPtr > m_alignment_sample_ptrs;
		uint32_t m_flanking_padding;
		uint32_t m_match_value;
		uint32_t m_mismatch_value;
		uint32_t m_gap_open_value;
		uint32_t m_gap_extension_value;
		ThreadPool m_thread_pool;
		bool m_print_graphs;
		int32_t m_mapping_quality;
		int32_t m_read_sample_limit;
		Sample::SharedPtr m_override_shared_ptr;
		std::mutex m_alignment_tracker_mutex;
		std::unordered_set< std::string > m_alignment_tracker_set;
	};
}

#endif //GRAPHITE_GRAPHPROCESSOR_H
