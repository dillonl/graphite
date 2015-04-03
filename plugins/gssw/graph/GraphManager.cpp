#include "core/variants/VariantList.h"
#include "core/utils/ThreadPool.hpp"
#include "GraphManager.h"
#include "AlignmentReporter.h"
#include "core/genotyper/IGenotyper.h"

#include "core/alignments/BamAlignmentReader.h"

#include <queue>


#include <boost/bind.hpp>

namespace gwiz
{
namespace gssw
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignmentReaderManager::SharedPtr alignmentReaderManager, IGraphAdjudicator::SharedPtr graphAdjudicatorPtr) :
		m_reference_ptr(referencePtr),
		m_variant_list_ptr(variantListPtr),
		m_alignment_reader_manager(alignmentReaderManager),
		m_graph_adjudicator_ptr(graphAdjudicatorPtr)
	{
	}

	IVariantList::SharedPtr GraphManager::buildGraphs(Region::SharedPtr regionPtr, size_t graphSize, size_t overlap, size_t alignmentPadding)
	{
		auto reportedVariants = std::make_shared< VariantList >();
		std::mutex reportedVariantsMutex;

		std::string referenceID = regionPtr->getReferenceID();
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();
		position currentPosition = startPosition;
		int64_t distance = endPosition - startPosition; // this number may become negative

		while (currentPosition < endPosition)
		{
			auto endGraphPosition = (currentPosition + graphSize > endPosition) ? endPosition : (currentPosition + graphSize);
			auto graphRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition) + "-" + std::to_string(endGraphPosition)));
			auto variantsListPtr = this->m_variant_list_ptr->getVariantsInRegion(graphRegion);
			if (variantsListPtr->getCount() == 0) { continue; }
			auto alignmentReaderPtr = this->m_alignment_reader_manager->generateAlignmentReader(); // create alignment reader

			auto alignmentRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition + alignmentPadding) + "-" + std::to_string(endGraphPosition - alignmentPadding)));
			// create region for alignmentReader
			alignmentReaderPtr->init();
			alignmentReaderPtr->setRegion(alignmentRegion); // set alignmentreader's region

			auto funct = std::bind(&GraphManager::constructAndAdjudicateGraph, this, reportedVariants, variantsListPtr, alignmentReaderPtr, std::ref(reportedVariantsMutex));
			ThreadPool::Instance()->enqueue(funct);
			currentPosition += graphSize - overlap;
		}
		return reportedVariants;
	}

	void GraphManager::constructAndAdjudicateGraph(IVariantList::SharedPtr reportedVariants, IVariantList::SharedPtr variantsListPtr, IAlignmentReader::SharedPtr alignmentReaderPtr, std::mutex& reportedVariantsMutex)
	{
		auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr);
		gsswGraph->constructGraph();
		auto variantsPtrList =  this->m_graph_adjudicator_ptr->adjudicateGraph(gsswGraph, alignmentReaderPtr);

		{
			std::lock_guard< std::mutex > lock(reportedVariantsMutex); // this is released when it falls out of scope
			reportedVariants->addVariants(variantsPtrList);
		}
	}


	/*
	void GraphManager::buildGraphs()
	{
	auto variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
		position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
		position previousVariantEndPosition = 0;
		Variant::SharedPtr variantPtr;
		Variant::SharedPtr firstVariantInList;
		size_t minContig = 500;
		size_t contig = 0;
		size_t maxContig = 0;
		position maxContigPosition = 0;
		size_t totalContigs = 0;
		std::map< size_t, size_t > contigCounterMap;
		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
			int32_t variantPositionDifference = (variantPtr->getPosition() - previousVariantEndPosition);
			if (firstVariantInList == NULL)
			{
				firstVariantInList = variantPtr;
			}
			uint32_t contigLength = (previousVariantEndPosition > 0) ? (previousVariantEndPosition - firstVariantInList->getPosition()) : this->m_padding;
			// if there is a gap between the variants that is larger than the padding then generate a graph and push it in the queue
			if ((contigLength > minContig) && (this->m_padding < variantPositionDifference))
			{
				++totalContigs;
				size_t contigCounter = 1;
				if (contigCounterMap.find(contigLength) != contigCounterMap.end())
				{
					contigCounter = contigCounterMap[contigLength] + 1;
				}
				contigCounterMap[contigLength] = contigCounter;
				if (maxContig < contigLength)
				{
					maxContig = contigLength;
					maxContigPosition = firstVariantInList->getPosition();
				}
				position startPosition = firstVariantInList->getPosition() - this->m_padding;
				position endPosition = previousVariantEndPosition;
				auto funct = std::bind(&GraphManager::buildGraph, this, startPosition, endPosition, variantListPtr);
				ThreadPool::Instance()->enqueue(funct);
				variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
				firstVariantInList = variantPtr; // set the first variant
				contigLength = 0;
			}
			variantListPtr->addVariant(variantPtr);
			previousVariantEndPosition = variantPtr->getPosition() + variantPtr->getRef()[0].size();
		}

		ThreadPool::Instance()->joinAll();
		AlignmentReporter::Instance()->printAlignmentReportsToStream(std::cout);
		// std::cout << "threads finished" << std::endl;
		// IGenotyper::Instance()->printVariants();

	}
*/

	IVariantList::SharedPtr GraphManager::buildGraph(position startPosition, position endPosition, IVariantList::SharedPtr variantListPtr)
	{
		/*
		auto reportedVariants = std::make_shared< VariantList >();
		static uint32_t contigsDone = 0;
		auto alignmentReaderPtr = this->m_alignment_reader_manager->generateAlignmentReader(); // create alignment reader

		// create region for alignmentReader
		auto region = std::make_shared< Region >(this->m_reference_ptr->getRegion()->getReferenceID() + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition));
		alignmentReaderPtr->init();
		alignmentReaderPtr->setRegion(region); // set alignmentreader's region
		auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantListPtr, alignmentReaderPtr);
		gsswGraph->constructGraph();
		return reportedVariants;
		*/
		return nullptr;
	}
}
}
