#include "GraphProcessor.h"
#include "ReferenceGraph.h"

#include <unordered_set>
#include <thread>
#include <queue>
#include <algorithm>
#include <functional>

namespace graphite
{
	GraphProcessor::GraphProcessor(FastaReference::SharedPtr fastaReferencePtr, const std::vector< BamReader::SharedPtr >& bamReaderPtrs, const std::vector< VCFReader::SharedPtr >& vcfReaderPtrs,  uint32_t matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue, bool printGraph, int32_t mappingQuality, int32_t readSampleLimit) :
		m_fasta_reference_ptr(fastaReferencePtr),
		m_bam_reader_ptrs(bamReaderPtrs),
		m_vcf_reader_ptrs(vcfReaderPtrs),
		m_flanking_padding(500),
		m_match_value(matchValue),
		m_mismatch_value(mismatchValue),
		m_gap_open_value(gapOpenValue),
		m_gap_extension_value(gapExtensionValue),
		m_thread_pool(std::thread::hardware_concurrency() * 2),
		m_print_graphs(printGraph),
		m_mapping_quality(mappingQuality),
		m_read_sample_limit(readSampleLimit),
		m_override_shared_ptr(nullptr)
	{
		for (auto bamReaderPtr : bamReaderPtrs)
		{
			auto bamSamplePtrs = bamReaderPtr->getSamplePtrs();
			for (auto samplePtr: bamSamplePtrs)
			{
				m_override_shared_ptr = samplePtr;
			}
		}
	}

	GraphProcessor::~GraphProcessor()
	{
	}

	void GraphProcessor::processVariants()
	{
		uint32_t graphSpacing = 200;
		bool overwriteSamples = false;
		// set the graph spacing to be the largest read size
		for (auto bamReaderPtr : this->m_bam_reader_ptrs)
		{
			graphSpacing = (graphSpacing >= bamReaderPtr->getReadLength()) ? graphSpacing : bamReaderPtr->getReadLength();
		}
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

		/*
		  graphPtr->printGraphVisOutput();
		auto graphPaths = graphPtr->getAllPathsAsStrings();
		std::cout << "[GraphProcessor.cpp] ------------"<< std::endl;
		for (auto path : graphPaths)
		{
			std::cout << path << std::endl;
		}
		std::cout << "[GraphProcessor]------------"<< std::endl;

		std::cout << "skipping adjudicating reads: GraphProcessor.cpp: 91" << std::endl;
		*/
		std::vector< Region::SharedPtr > graphRegionPtrs = graphPtr->getRegionPtrs();

		// get all alignments
		std::vector< std::shared_ptr< BamAlignment > > bamAlignmentPtrs;

		getAlignmentsInRegion(bamAlignmentPtrs, graphRegionPtrs, true);
		for (auto bamAlignmentPtr : bamAlignmentPtrs)
		{
			std::string sampleName;
			Sample::SharedPtr samplePtr = m_override_shared_ptr;
			if (m_override_shared_ptr == nullptr)
			{
				bamAlignmentPtr->GetTag("RG", sampleName);
			}
			else
			{
				sampleName = samplePtr->getName();
			}
			auto iter = this->m_bam_sample_ptrs.find(sampleName);
			if (iter != this->m_bam_sample_ptrs.end() || samplePtr != nullptr)
			{
				if (samplePtr == nullptr)
				{
					samplePtr = iter->second;
				}
				uint32_t matchValue = m_match_value;
				uint32_t mismatchValue = m_mismatch_value;
				uint32_t gapOpenValue = m_gap_open_value;
				uint32_t gapExtensionValue = m_gap_extension_value;
				auto funct = [graphPtr,bamAlignmentPtr, samplePtr, matchValue, mismatchValue, gapOpenValue, gapExtensionValue]()
				{
					graphPtr->adjudicateAlignment(bamAlignmentPtr, samplePtr, matchValue, mismatchValue, gapOpenValue, gapExtensionValue, 0);
				};
				m_thread_pool.enqueue(funct);
			}
		}
		m_thread_pool.join();
		bamAlignmentPtrs.clear();
		graphPtr->clearResources();
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
					bamReaderPtr->fetchBamAlignmentPtrsInRegion(bamAlignmentPtrsTmp, flankingRegionPtr, false, false, this->m_mapping_quality);
				}
				bamReaderPtr->fetchBamAlignmentPtrsInRegion(bamAlignmentPtrsTmp, regionPtr, false, false, this->m_mapping_quality);
				if (iter == regionPtrs.end() && getFlankingUnalignedReads)
				{
					auto flankingRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), regionPtr->getEndPosition(), regionPtr->getEndPosition() + this->m_flanking_padding, regionPtr->getBased());
					bamReaderPtr->fetchBamAlignmentPtrsInRegion(bamAlignmentPtrsTmp, flankingRegionPtr, false, false, this->m_mapping_quality);
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
		if (this->m_read_sample_limit > 0)
		{
			std::shuffle(bamAlignmentPtrs.begin(), bamAlignmentPtrs.end(), default_random_engine(0));
			bamAlignmentPtrs.resize(this->m_read_sample_limit);
		}
	}
}
