#include "GraphManager.h"
#include "ReferenceGraph.h"
#include "core/alignment/AlignmentReporter.h"
#include "core/util/ThreadPool.hpp"
#include "core/variant/VariantList.h"
#include "core/mapping/MappingManager.h"
#include "core/mapping/GSSWMapping.h"

#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignment.h"
//#include "core/file/FastaFileWriter.h"

#include <queue>
#include <algorithm>
#include <functional>

//#include <string>

namespace graphite
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantManager::SharedPtr variantManagerPtr, IAlignmentManager::SharedPtr alignmentManagerPtr, IAdjudicator::SharedPtr adjudicatorPtr) :
		m_reference_ptr(referencePtr),
		m_variant_manager_ptr(variantManagerPtr),
		m_alignment_manager_ptr(alignmentManagerPtr),
		m_adjudicator_ptr(adjudicatorPtr)
	{
	}

	void GraphManager::buildGraphs(Region::SharedPtr regionPtr, uint32_t readLength)
	{
		auto variantsListPtr = this->m_variant_manager_ptr->getVariantsInRegion(regionPtr);
		if (variantsListPtr->getCount() == 0) // if we don't have variants or alignments in the region, then return
		{
			return;
        }

		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		// loop through variants and build and adjudicate graphs
		IVariant::SharedPtr variantPtr = nullptr;
		while (variantsListPtr->getNextVariant(variantPtr))
		{
			position startPosition = 0;
			position endPosition = 0;
			auto variantRegionPtrs = variantPtr->getRegions();
			if (variantRegionPtrs.size() > 0)
			{
				startPosition = variantPtr->getRegions()[0]->getStartPosition();
				endPosition = variantPtr->getRegions()[variantPtr->getRegions().size() - 1]->getEndPosition();
			}
			std::vector< IVariant::SharedPtr > variantPtrs;
			IVariant::SharedPtr nextVariantPtr = nullptr;
			variantPtrs.emplace_back(variantPtr);
			while (variantsListPtr->peekNextVariant(nextVariantPtr) && variantPtr->doesOverlap(nextVariantPtr) && (!variantPtr->isStructuralVariant() && !nextVariantPtr->isStructuralVariant()))
			{
				variantsListPtr->getNextVariant(nextVariantPtr);
				variantPtrs.emplace_back(nextVariantPtr);

				for (auto regionPtr : nextVariantPtr->getRegions())
				{
					auto tmpStartPosition = regionPtr->getStartPosition();
					auto tmpEndPosition = regionPtr->getEndPosition();
					startPosition = (tmpStartPosition < startPosition) ? tmpStartPosition : startPosition;
					endPosition = (tmpEndPosition > endPosition) ? tmpEndPosition : endPosition;
				}
			}
			// getting all the alignments in the variant's region
			std::vector< IAlignment::SharedPtr > alignmentPtrs;
			std::unordered_set< IAlignment* > alignmentPtrSet;
			for (auto regionVariantPtr : variantPtrs)
			{
				auto regionPtrs = regionVariantPtr->getRegions();
				for (auto regionPtr : regionPtrs)
				{
					auto regionAlignmentListPtr = this->m_alignment_manager_ptr->getAlignmentsInRegion(regionPtr);
					auto regionAlignmentPtrsCount = regionAlignmentListPtr->getCount();
					if (regionAlignmentPtrsCount > 0)
					{
						auto regionAlignmentPtrs = regionAlignmentListPtr->getAlignmentPtrs();
						// alignmentPtrs.insert(alignmentPtrs.end(), regionAlignmentPtrs.begin(), regionAlignmentPtrs.end()); // combine lists of alignments
						for (auto alignmentPtr : regionAlignmentPtrs)
						{
							if (alignmentPtrSet.find(alignmentPtr.get()) == alignmentPtrSet.end())
							{
								alignmentPtrSet.emplace(alignmentPtr.get());
								alignmentPtrs.emplace_back(alignmentPtr);
								if (alignmentPtr->getPosition() < startPosition) { startPosition = alignmentPtr->getPosition(); }
								if ((alignmentPtr->getPosition() + alignmentPtr->getLength()) > endPosition) { endPosition = alignmentPtr->getPosition() + alignmentPtr->getLength(); }
							}
						}
					}
				}
			}
			alignmentPtrSet.clear();
			if (alignmentPtrs.size() > 0)
			{
				// find the start and end position for the graph
				auto graphAlignmentRegion = std::make_shared< Region >(regionPtr->getReferenceID(), startPosition, endPosition, Region::BASED::ONE);
				for (auto i = 1; i < alignmentPtrs.size(); ++i)
				{
					if (alignmentPtrs[i]->getPosition() < startPosition) { startPosition = alignmentPtrs[i]->getPosition(); }
					if (alignmentPtrs[i]->getPosition() + alignmentPtrs[0]->getLength() > endPosition) { endPosition = alignmentPtrs[i]->getPosition() + alignmentPtrs[0]->getLength(); }
				}
				auto variantListPtr = std::make_shared< VariantList >(variantPtrs, this->m_reference_ptr);
				auto alignmentListPtr = std::make_shared< AlignmentList >(alignmentPtrs);
				constructAndAdjudicateGraph(variantListPtr, alignmentListPtr, graphAlignmentRegion, readLength);
			}
		}

		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) != std::future_status::ready)
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}
	}

    void GraphManager::registerNodeInfo(uint32_t nodeID, int32_t length, NodeInfo::VariantType variantType)
    {
        auto nodeInfoPtr = std::make_shared< NodeInfo >(length, variantType);
        m_node_info_map.insert( {nodeID, nodeInfoPtr} );
    }

	void GraphManager::constructAndAdjudicateGraph(IVariantList::SharedPtr variantsListPtr, IAlignmentList::SharedPtr alignmentListPtr, Region::SharedPtr regionPtr, uint32_t readLength)
	{
		uint32_t numGraphCopies = (alignmentListPtr->getCount() < ThreadPool::Instance()->getThreadCount()) ? alignmentListPtr->getCount() : ThreadPool::Instance()->getThreadCount();  // get the min of threadcount and alignment count, this is the num of simultanious threads processing this graph
		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;

		auto gsswGraphPtr = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, regionPtr, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue(), numGraphCopies);
		gsswGraphPtr->constructGraph();

		auto referenceGraphPtr = std::make_shared< ReferenceGraph >(this->m_reference_ptr, variantsListPtr, regionPtr, this->m_adjudicator_ptr->getMatchValue(), this->m_adjudicator_ptr->getMisMatchValue(), this->m_adjudicator_ptr->getGapOpenValue(), this->m_adjudicator_ptr->getGapExtensionValue(), numGraphCopies);
		referenceGraphPtr->constructGraph();

		static int count = 0;

		IAlignment::SharedPtr alignmentPtr;
		auto alignmentPtrs = alignmentListPtr->getAlignmentPtrs();

		while (alignmentListPtr->getNextAlignment(alignmentPtr))
		{
            position variantPosition = gsswGraphPtr->getVariantPosition();

			auto gsswGraphContainer = gsswGraphPtr->getGraphContainer();
			auto refGraphContainer = referenceGraphPtr->getGraphContainer();
			auto funct = [gsswGraphContainer, refGraphContainer, gsswGraphPtr, referenceGraphPtr, alignmentPtr, this, variantPosition]()
			{
                std::string originalGraphPathHeader;
                std::string alteredGraphPathHeader;
                std::string originalGraphPathSequence;
                std::string alteredGraphPathSequence;
                
				auto refTraceback = referenceGraphPtr->traceBackAlignment(alignmentPtr, refGraphContainer);
				auto referenceMappingPtr = std::make_shared< GSSWMapping >(refTraceback, alignmentPtr);
				auto referenceSWScore = referenceMappingPtr->getMappingScore();
				uint32_t referenceSWPercent = ((referenceSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);

				auto tracebackPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr, gsswGraphContainer);
				auto gsswMappingPtr = std::make_shared< GSSWMapping >(tracebackPtr, alignmentPtr);

                // Write Bam data to new sam file.
                // If I can tie the header to the appropriate alignment then I can use the same process to tie the new position to the appropriate alignment.
                // Verify:
                //   New position caclulation (col 3)
                //   New CIGAR string (col 7)
                // Remember that SAM has to be sorted by position to be loaded into IGV.
                {
                    std::lock_guard< std::mutex > lock(this->m_gssw_graph_mutex);
                    // This function also updates the bamAlignment with the proper read group.
                    referenceMappingPtr->getGraphPathHeaderAndSequence(originalGraphPathHeader, originalGraphPathSequence, variantPosition);
                    gsswMappingPtr->getGraphPathHeaderAndSequence(alteredGraphPathHeader, alteredGraphPathSequence, variantPosition);
                    m_header_sequence_map.insert( {originalGraphPathHeader, originalGraphPathSequence} ); 
                    m_header_sequence_map.insert( {alteredGraphPathHeader, alteredGraphPathSequence} );

                    // Get data for NodeInfo
                    std::vector< uint32_t > originalNodeIDs = referenceMappingPtr->getNodeIDs();
                    std::vector< int32_t > originalNodeLengths = referenceMappingPtr->getNodeLengths();
                    std::vector< NodeInfo::VariantType > originalVariantTypes = gsswMappingPtr->getVariantTypes();
                    std::vector< uint32_t > alteredNodeIDs = gsswMappingPtr->getNodeIDs();
                    std::vector< int32_t > alteredNodeLengths = gsswMappingPtr->getNodeLengths();
                    std::vector< NodeInfo::VariantType > alteredVariantTypes = gsswMappingPtr->getVariantTypes();
                    
                    for (int i = 0; i < originalNodeIDs.size(); ++i)
                    {
                        registerNodeInfo(originalNodeIDs[i], originalNodeLengths[i], originalVariantTypes[i]);
                        registerNodeInfo(alteredNodeIDs[i], alteredNodeLengths[i], alteredVariantTypes[i]);
                    }

                    graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(alignmentPtr);
                    std::string readGroup;
                    if (gsswMappingPtr->getAltCount() > 0)
                        readGroup = "REF";
                    else
                        readGroup = "ALT"; 
                    //m_bam_alignment_ptrs.push_back(alignmentPtr);
                    //graphite::BamAlignment::SharedPtr alteredBamAlignmentPtr;

                    std::ofstream samFile;
                    samFile.open("TempAlignmentFile.sam", std::ios::app);

                    // Write out original bamAlignment.
                    samFile
                        << bamAlignmentPtr->getName()                   << "\t" //  1. QNAME
                        << bamAlignmentPtr->getAlignmentFlag()          << "\t" //  2. FLAG
                        << originalGraphPathHeader                      << "\t" //  3. RNAME
                        // Need to find out why I need to + 1 on the offset.
                        << referenceMappingPtr->getOffset()             << "\t" //  4. POS New position.
                        << bamAlignmentPtr->getOriginalMapQuality()     << "\t" //  5. MAPQ
                        << bamAlignmentPtr->getCigarString()            << "\t" //  6. New CIGAR string.
                        << bamAlignmentPtr->getMateReferenceName()      << "\t" //  7. Place holder for actual value.
                        << bamAlignmentPtr->getMatePosition() + 1       << "\t" //  8. PNEXT +1 because BamTools mate position is 0-based.
                        << bamAlignmentPtr->getTemplateLength()         << "\t" //  9. TLEN
                        << bamAlignmentPtr->getSequence()               << "\t" // 10. SEQ
                        << bamAlignmentPtr->getFastqQualities()         << "\t" // 11. QUAL
                        << "RG:Z:" << readGroup                                 // 12. Optional read group field.
                        << std::endl; 

                    // Write out updated bamAlignment.
                    samFile
                        << bamAlignmentPtr->getName()                   << "\t" //  1. QNAME
                        << bamAlignmentPtr->getAlignmentFlag()          << "\t" //  2. FLAG
                        << alteredGraphPathHeader                       << "\t" //  3. RNAME;
                        // Need to find out why I need to + 1 on the offset.
                        << gsswMappingPtr->getOffset()                  << "\t" //  4. POS New position.
                        << bamAlignmentPtr->getOriginalMapQuality()     << "\t" //  5. MAPQ
                        << gsswMappingPtr->getCigarString(m_adjudicator_ptr) << "\t" //  6. New CIGAR string.
                        << bamAlignmentPtr->getMateReferenceName()      << "\t" //  7. RNEXT INCORRECT value.
                        << bamAlignmentPtr->getMatePosition() + 1       << "\t" //  8. PNEXT +1 because BamTools mate position is 0-based.
                        << bamAlignmentPtr->getTemplateLength()         << "\t" //  9. TLEN
                        << bamAlignmentPtr->getSequence()               << "\t" // 10. SEQ
                        << bamAlignmentPtr->getFastqQualities()         << "\t" // 11. QUAL
                        << "RG:Z:" << readGroup                                 // 12. Optional read group field.
                        << std::endl; 

                    samFile.close();     // Close the SAM file.
                }

				auto gsswSWScore = referenceMappingPtr->getMappingScore();
				uint32_t gsswSWPercent = ((gsswSWScore / (double)(alignmentPtr->getLength() * this->m_adjudicator_ptr->getMatchValue())) * 100);

				if (this->m_adjudicator_ptr->adjudicateMapping(gsswMappingPtr, referenceSWPercent))
				{
					MappingManager::Instance()->registerMapping(gsswMappingPtr);
				}
		    };

			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
		}

		// wait for all the functions to complete (for all graphs to be created and adjudicated)
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) != std::future_status::ready)
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}
	}

    std::unordered_map< std::string, std::string > GraphManager::getHeaderSequenceMap () { return m_header_sequence_map; }
    std::unordered_map< uint32_t, NodeInfo::SharedPtr > GraphManager::getNodeInfoMap () { return m_node_info_map; }
}
