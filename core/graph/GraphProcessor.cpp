#include "GraphProcessor.h"
#include "ReferenceGraph.h"
#include "GraphTraceback.hpp"

#include <unordered_set>
#include <thread>
#include <queue>
#include <algorithm>
#include <functional>
#include <random>
#include <atomic>

namespace graphite
{
	GraphProcessor::GraphProcessor(FastaReference::SharedPtr fastaReferencePtr, const std::vector< AlignmentReader::SharedPtr >& alignmentReaderPtrs, const std::vector< VCFReader::SharedPtr >& vcfReaderPtrs,  uint32_t matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue, bool printGraph, int32_t mappingQuality, int32_t readSampleLimit, uint32_t numberOfThreads) :
		m_fasta_reference_ptr(fastaReferencePtr),
		m_alignment_reader_ptrs(alignmentReaderPtrs),
		m_vcf_reader_ptrs(vcfReaderPtrs),
		m_flanking_padding(1),
		m_match_value(matchValue),
		m_mismatch_value(mismatchValue),
		m_gap_open_value(gapOpenValue),
		m_gap_extension_value(gapExtensionValue),
		m_thread_pool(numberOfThreads),
		m_print_graphs(printGraph),
		m_mapping_quality(mappingQuality),
		m_read_sample_limit(readSampleLimit),
		m_override_shared_ptr(nullptr)
	{
	}

	GraphProcessor::~GraphProcessor()
	{
	}

	void GraphProcessor::processVariants()
	{
		uint32_t graphSpacing = 10;
		bool overwriteSamples = false;
		// set the graph spacing to be the largest read size
		for (auto alignmentReaderPtr : this->m_alignment_reader_ptrs)
		{
			graphSpacing = (graphSpacing >= alignmentReaderPtr->getReadLength()) ? graphSpacing : alignmentReaderPtr->getReadLength();
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
		std::vector< Region::SharedPtr > graphRegionPtrs = graphPtr->getRegionPtrs();

		// get all alignments
		std::vector< Alignment::SharedPtr > alignmentPtrs;

		getAlignmentsInRegion(alignmentPtrs, graphRegionPtrs, true);
		for (auto alignmentPtr : alignmentPtrs)
		{
			std::string sampleName;
			sampleName = alignmentPtr->getSample()->getName();
			uint32_t matchValue = m_match_value;
			uint32_t mismatchValue = m_mismatch_value;
			uint32_t gapOpenValue = m_gap_open_value;
			uint32_t gapExtensionValue = m_gap_extension_value;
			// check if the alignment has already been processed
			/*
			{
				std::lock_guard< std::mutex > l(m_alignment_tracker_mutex);
				if (m_alignment_tracker_set.find(alignmentPtr->getUniqueReadName()) != m_alignment_tracker_set.end())
				{
					continue;
				}
				m_alignment_tracker_set.emplace(alignmentPtr->getUniqueReadName());
			}
			*/
			auto funct = [graphPtr, alignmentPtr, matchValue, mismatchValue, gapOpenValue, gapExtensionValue, variantPtrs]()
				{
					auto graphTraceback = std::make_shared< GraphTraceback >(graphPtr, matchValue, mismatchValue, gapOpenValue, gapExtensionValue);
					graphTraceback->processGraph(alignmentPtr);
					auto tracebackNodePtrs = graphTraceback->getTracebackNodePtrs();
					auto originalGraphTracebackCigar = graphTraceback->getNormalizedCigarString();
					int origTracebackSoftclipCount = std::count(originalGraphTracebackCigar.begin(), originalGraphTracebackCigar.end(), 'S');
					int counter = -1;
					if (graphTraceback->getTotalScore() >= 90)
					{
						for (auto nodePtr : tracebackNodePtrs)
						{

							if (!nodePtr->hasSiblings()) // if this is a reference "backbone" node that connects variants then skip it
							{
								continue;
							}

							auto altGraphPtr = graphPtr->createCopy();
							altGraphPtr->removeNodePtr(nodePtr);
							auto altGraphTraceback = std::make_shared< GraphTraceback >(altGraphPtr, matchValue, mismatchValue, gapOpenValue, gapExtensionValue);
							altGraphTraceback->processGraph(alignmentPtr);
							auto altGraphTracebackCigar = altGraphTraceback->getNormalizedCigarString();
							auto nodeScorePercent = graphTraceback->getNodeScorePercent(nodePtr);
							auto cigarComparisonEqual = originalGraphTracebackCigar.compare(altGraphTracebackCigar);
							int altTracebackSoftclipCount = std::count(altGraphTracebackCigar.begin(), altGraphTracebackCigar.end(), 'S');

							// static std::mutex l;
							// std::lock_guard< std::mutex > lock(l);
							// std::cout << "---------------" << std::endl;
							// std::cout << "orig: " << graphTraceback->getTracebackAsSequence("|") << "\tscore: " << graphTraceback->getTotalScore() << std::endl;
							// std::cout << "alt:  " << altGraphTraceback->getTracebackAsSequence("|") << "\tscore: " << altGraphTraceback->getTotalScore() << std::endl;
							// std::cout << "\torigCig: " << originalGraphTracebackCigar << std::endl;
							// std::cout << "\taltCig:  " << altGraphTracebackCigar << std::endl;
							// if (cigarComparisonEqual == 0)
							// {
								// std::cout << "PASSED - CIGAR MATCH" << std::endl;
							// }
							// if (graphTraceback->getTotalScore() == altGraphTraceback->getTotalScore())
							// {
								// std::cout << "scores are equal" << std::endl;
							// }

							if (nodeScorePercent >= 70)
							{
								for (auto nodeAllelePtr : nodePtr->getAllelePtrs())
								{
									//if (graphTraceback->getTotalScore() == altGraphTraceback->getTotalScore()) // check the cigar here if you want to
									if ((cigarComparisonEqual == 0 || graphTraceback->getTotalScore() == altGraphTraceback->getTotalScore()) && (origTracebackSoftclipCount == altTracebackSoftclipCount))
									{
										// this is if the node with that alignment is ambiguous
										nodeAllelePtr->incrementScoreCount(alignmentPtr, -1);
									}
									else
									{
										/*
										if (nodePtr->getAlleleType() == Node::ALLELE_TYPE::REF)
										{
											std::cout << "---------------" << std::endl;

											std::cout << "all paths: " << std::endl;
											auto paths = altGraphPtr->getAllPathsAsStrings();
											for (auto path : paths)
											{
												std::cout << "path: " << path << std::endl;
											}

											std::cout << "alignment: " << alignmentPtr->getSequence() << std::endl;
											std::cout << "we're counting allele: " << nodeAllelePtr->getSequence() << std::endl;
											std::cout << "orig: " << graphTraceback->getTracebackAsSequence("|") << "\tscore: " << graphTraceback->getTotalScore() << std::endl;
											std::cout << "alt:  " << altGraphTraceback->getTracebackAsSequence("|") << "\tscore: " << altGraphTraceback->getTotalScore() << std::endl;
											std::cout << "\torigCig: " << originalGraphTracebackCigar << std::endl;
											std::cout << "\taltCig:  " << altGraphTracebackCigar << std::endl;
											std::cout << "---------------" << std::endl;
										}
										*/
										nodeAllelePtr->incrementScoreCount(alignmentPtr, nodeScorePercent);
									}
								}
							}
							// std::cout << "---------------" << std::endl;
						}
					}
				};
			m_thread_pool.enqueue(funct);
		}
		m_thread_pool.join();
	}

	void GraphProcessor::getAlignmentsInRegion(std::vector< Alignment::SharedPtr >& alignmentPtrs, std::vector< Region::SharedPtr > regionPtrs, bool getFlankingUnalignedReads)
	{
		std::vector< Alignment::SharedPtr > alignmentPtrsTmp;
		for (auto iter = regionPtrs.begin(); iter != regionPtrs.end(); ++iter)
		{
			auto regionPtr = (*iter);
			for (auto alignmentReaderPtr : this->m_alignment_reader_ptrs)
			{
				if (iter == regionPtrs.begin() && getFlankingUnalignedReads)
				{
					auto flankingRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), regionPtr->getStartPosition() - this->m_flanking_padding, regionPtr->getStartPosition(), regionPtr->getBased());
					alignmentReaderPtr->fetchAlignmentPtrsInRegion(alignmentPtrsTmp, flankingRegionPtr, false, false, this->m_mapping_quality);
				}
				alignmentReaderPtr->fetchAlignmentPtrsInRegion(alignmentPtrsTmp, regionPtr, false, false, this->m_mapping_quality);
				if (iter == regionPtrs.end() && getFlankingUnalignedReads)
				{
					auto flankingRegionPtr = std::make_shared< Region >(regionPtr->getReferenceID(), regionPtr->getEndPosition(), regionPtr->getEndPosition() + this->m_flanking_padding, regionPtr->getBased());
					alignmentReaderPtr->fetchAlignmentPtrsInRegion(alignmentPtrsTmp, flankingRegionPtr, false, false, this->m_mapping_quality);
				}
			}
		}

		// make sure the reads only appear in alignmentPtrs once
		alignmentPtrs.clear();
		std::unordered_set< std::string > alignmentIDTracker;
		for (auto alignmentPtr : alignmentPtrsTmp)
		{
			std::string alignmentID = alignmentPtr->getUniqueReadName();
			if (alignmentIDTracker.find(alignmentID) == alignmentIDTracker.end())
			{
				alignmentPtrs.emplace_back(alignmentPtr);
				alignmentIDTracker.emplace(alignmentID);
			}
		}
		if (this->m_read_sample_limit < alignmentPtrs.size())
		{
			std::shuffle(alignmentPtrs.begin(), alignmentPtrs.end(), default_random_engine(0));
			alignmentPtrs.resize(this->m_read_sample_limit);
		}
	}
}
