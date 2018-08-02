#include "GraphProcessor.h"

#include <unordered_set>
#include <thread>
#include <queue>
#include <algorithm>
#include <functional>

namespace graphite
{
	GraphProcessor::GraphProcessor(FastaReference::SharedPtr fastaReferencePtr, const std::vector< BamReader::SharedPtr >& bamReaderPtrs, const std::vector< VCFReader::SharedPtr >& vcfReaderPtrs,  uint32_t matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue, bool printGraph) :
		m_fasta_reference_ptr(fastaReferencePtr),
		m_bam_reader_ptrs(bamReaderPtrs),
		m_vcf_reader_ptrs(vcfReaderPtrs),
		m_flanking_padding(500),
		m_match_value(matchValue),
		m_mismatch_value(mismatchValue),
		m_gap_open_value(gapOpenValue),
		m_gap_extension_value(gapExtensionValue),
		m_thread_pool(std::thread::hardware_concurrency()),
		m_print_graphs(printGraph)
	{
		for (auto bamReaderPtr : bamReaderPtrs)
		{
			auto bamSamplePtrs = bamReaderPtr->getSamplePtrs();
			for (auto samplePtr: bamSamplePtrs)
			{
				this->m_bam_sample_ptrs.emplace(samplePtr->getReadgroup(), samplePtr);
			}
		}
	}

	GraphProcessor::~GraphProcessor()
	{
	}

	void GraphProcessor::processVariants()
	{
		uint32_t graphSpacing = 200;
		// set the graph spacing to be the largest read size
		for (auto bamReaderPtr : this->m_bam_reader_ptrs) { graphSpacing = (graphSpacing >= bamReaderPtr->getReadLength()) ? graphSpacing : bamReaderPtr->getReadLength(); }
		std::vector< Variant::SharedPtr > variantPtrs;
		while (true)
		{
			variantPtrs.clear();
			for (auto vcfReaderPtr : this->m_vcf_reader_ptrs)
			{
				vcfReaderPtr->getNextVariants(variantPtrs, graphSpacing);
			}
			// only adjudicate if there are variants to adjudicate
			if (variantPtrs.size() > 0)
			{
				adjudicateVariants(variantPtrs, graphSpacing);
				// for (auto variantPtr : variantPtrs)
				for (int i = 0; i < variantPtrs.size(); ++i)
				{
					auto variantPtr = variantPtrs[i];
					variantPtr->writeVariant();
				}
			}
			else
			{
				break;
			}
		}
	}

	void GraphProcessor::adjudicateVariants(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t graphSpacing)
	{
		// generate graph
		auto graphPtr = std::make_shared< Graph >(this->m_fasta_reference_ptr, variantPtrs, graphSpacing, this->m_print_graphs);
		std::vector< Region::SharedPtr > graphRegionPtrs = graphPtr->getRegionPtrs();

		// get all alignments
		std::vector< std::shared_ptr< BamAlignment > > bamAlignmentPtrs;

		getAlignmentsInRegion(bamAlignmentPtrs, graphRegionPtrs, true);
		for (auto bamAlignmentPtr : bamAlignmentPtrs)
		{
			std::string sampleName;
			bamAlignmentPtr->GetTag("RG", sampleName);
			auto iter = this->m_bam_sample_ptrs.find(sampleName);
			if (iter != this->m_bam_sample_ptrs.end())
			{
				// graphPtr->adjudicateAlignment(bamAlignmentPtr, iter->second, m_match_value, m_mismatch_value, m_gap_open_value, m_gap_extension_value);
				auto funct = std::bind(&Graph::adjudicateAlignment, graphPtr, bamAlignmentPtr, iter->second, m_match_value, m_mismatch_value, m_gap_open_value, m_gap_extension_value);
				m_thread_pool.enqueue(funct);
			}
		}
		m_thread_pool.join();
		bamAlignmentPtrs.clear();
	}

	void GraphProcessor::adjudicateVariants2(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t graphSpacing)
	{
		// generate graph
		auto graphPtr = std::make_shared< Graph >(this->m_fasta_reference_ptr, variantPtrs, graphSpacing, this->m_print_graphs);
		std::vector< Region::SharedPtr > graphRegionPtrs = graphPtr->getRegionPtrs();

		// get all alignments
		std::vector< std::shared_ptr< BamAlignment > > bamAlignmentPtrs;
		getAlignmentsInRegion(bamAlignmentPtrs, graphRegionPtrs, true);

        // uint32_t numThreads = std::thread::hardware_concurrency() * 2;
		/*
		std::cout << "-----start-----" << std::endl;
		std::cout << "adjudicating: " << bamAlignmentPtrs.size() << " alignments" << std::endl;
		for (auto regionPtr : graphRegionPtrs)
		{
			std::cout << "region: " << regionPtr->getStartPosition() << " " << regionPtr->getEndPosition() << std::endl;
		}
		*/
		std::vector< std::future< void > > futures;
		// uint32_t numThreads = 2;
		for (auto bamAlignmentPtr : bamAlignmentPtrs)
		{

			std::string sampleName;
			bamAlignmentPtr->GetTag("RG", sampleName);
			auto iter = this->m_bam_sample_ptrs.find(sampleName);
			if (iter != this->m_bam_sample_ptrs.end())
			{
				// futures.emplace_back(this->m_threadpool.enqueue(std::bind(&Graph::adjudicateAlignment, graphPtr, bamAlignmentPtr, iter->second, m_match_value, m_mismatch_value, m_gap_open_value, m_gap_extension_value)));
			}
		}
		for (auto& f : futures)
		{
			f.wait();
		}
		futures.clear();
		bamAlignmentPtrs.clear();
		// std::cout << "-----end-----" << std::endl;
	}

	void GraphProcessor::getAlignmentsInRegion(std::vector< std::shared_ptr< BamAlignment > >& bamAlignmentPtrs, std::vector< Region::SharedPtr > regionPtrs, bool getFlankingUnalignedReads)
	{
		std::vector< std::shared_ptr< BamAlignment > > bamAlignmentPtrsTmp;
		for (auto iter = regionPtrs.begin(); iter != regionPtrs.end(); ++iter)
		{
			auto regionPtr = (*iter);
			for (auto bamReaderPtr : this->m_bam_reader_ptrs)
			{
				if (iter == regionPtrs.begin() && getFlankingUnalignedReads)
				{
					auto flankingRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), regionPtr->getStartPosition() - this->m_flanking_padding, regionPtr->getStartPosition(), regionPtr->getBased());
					bamReaderPtr->fetchBamAlignmentPtrsInRegion(bamAlignmentPtrsTmp, flankingRegionPtr, false, false);
				}
				bamReaderPtr->fetchBamAlignmentPtrsInRegion(bamAlignmentPtrsTmp, regionPtr, false, false);
				if (iter == regionPtrs.end() && getFlankingUnalignedReads)
				{
					auto flankingRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), regionPtr->getEndPosition(), regionPtr->getEndPosition() + this->m_flanking_padding, regionPtr->getBased());
					bamReaderPtr->fetchBamAlignmentPtrsInRegion(bamAlignmentPtrsTmp, flankingRegionPtr, false, false);
				}
			}
		}

		// make sure the reads only appear in bamAlignmentPtrs once
		bamAlignmentPtrs.clear();
		std::unordered_set< std::string > bamIDTracker;
		for (auto bamAlignmentPtr : bamAlignmentPtrsTmp)
		{
			std::string bamID = bamAlignmentPtr->Name + std::to_string(bamAlignmentPtr->IsFirstMate());
			if (bamIDTracker.find(bamID) == bamIDTracker.end())
			{
				bamAlignmentPtrs.emplace_back(bamAlignmentPtr);
				bamIDTracker.emplace(bamID);
			}
		}
	}
}
