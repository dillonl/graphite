#include "core/variants/VariantList.h"
#include "core/utils/ThreadPool.hpp"
#include "GraphManager.h"
#include "core/genotyper/IGenotyper.h"

#include "core/alignments/BamAlignmentReader.h"

#include <queue>


#include <boost/bind.hpp>

namespace gwiz
{
namespace gssw
{

	GraphManager::GraphManager(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignmentReaderManager::SharedPtr alignmentReaderManager, size_t padding) :
		m_reference_ptr(referencePtr),
		m_variant_list_ptr(variantListPtr),
		m_alignment_reader_manager(alignmentReaderManager),
		m_padding(padding)
	{
	}

	void GraphManager::buildGraphs()
	{
		auto variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
		position previousVariantEndPosition = 0;
		Variant::SharedPtr variantPtr;
		Variant::SharedPtr firstVariantInList;
		size_t minContig = 500;
		size_t contig = 0;
		size_t count = 0;
		size_t maxContig = 0;
		position maxContigPosition = 0;
		size_t totalContigs = 0;
		std::map< size_t, size_t > contigCounterMap;
		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
			++count;
			int32_t variantPositionDifference = (variantPtr->getPosition() - previousVariantEndPosition);
			if (firstVariantInList == NULL)
			{
				firstVariantInList = variantPtr;
			}
			uint32_t contigLength = (previousVariantEndPosition > 0) ? (previousVariantEndPosition - firstVariantInList->getPosition()) : this->m_padding;
			/*
			if (count % 10000 == 0)
			{
				std::cout << count << std::endl;
			}
			*/
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
		/*
		for (auto value : contigCounterMap)
		{
			std::cout << "Contig Size: " << std::get< 0 >(value) << " " << std::get< 1 >(value) << std::endl;
		}
		*/

		std::cout << "Max Contig: " << maxContig << std::endl;
		std::cout << "Max Contig Position: " << maxContigPosition << std::endl;
		std::cout << "total contigs: " << totalContigs << std::endl;
		std::cout << "finished" << std::endl;
		ThreadPool::Instance()->joinAll();
		std::cout << "threads finished" << std::endl;
		IGenotyper::Instance()->printVariants();

	}

	void GraphManager::buildGraph(position startPosition, position endPosition, IVariantList::SharedPtr variantListPtr)
	{
		// std::cout << startPosition << " " << endPosition << std::endl;
		// int x = 0;
		// std::cin >> x;
		static uint32_t contigsDone = 0;
		auto alignmentReaderPtr = this->m_alignment_reader_manager->generateAlignmentReader(); // create alignment reader

		// create region for alignmentReader
		auto region = std::make_shared< Region >(this->m_reference_ptr->getRegion()->getReferenceID() + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition));
		alignmentReaderPtr->init();
		alignmentReaderPtr->setRegion(region); // set alignmentreader's region
		auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantListPtr, alignmentReaderPtr);
		gsswGraph->constructGraph();
		// this->m_gssw_graphs.push(gsswGraph);
		// if (region->getEndPosition() == 36927311)
		// if ( contigsDone >= 97310)
		// {
		// 	GSSWGraph::PrintStuff = true;
		// }
		/*
		if (++contigsDone % 1000 == 0)
		{
			std::cout << "Contigs Computed: " << contigsDone << std::endl;
			// auto region = alignmentReaderPtr->getRegion();
			std::cout << "region: " << region->getStartPosition() << " " << region->getEndPosition() << std::endl;
		}
		*/
	}
}
}
