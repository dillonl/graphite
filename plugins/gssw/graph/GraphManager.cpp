#include "core/variants/VariantList.h"
#include "core/utils/ThreadPool.hpp"
#include "GraphManager.h"

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
		size_t minContig = 1000;
		size_t contig = 0;
		size_t count = 0;
		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
			++count;
			int32_t variantPositionDifference = (variantPtr->getPosition() - previousVariantEndPosition);
			if (firstVariantInList == NULL)
			{
				firstVariantInList = variantPtr;
			}
			uint32_t contigLength = (previousVariantEndPosition - firstVariantInList->getPosition());
			// if (count >= 86458)
			// if (BamAlignmentReader::ReaderCount > 10)
			// while (BamAlignmentReader::ReaderCount > 100)
			// {
				// std::cout <<  BamAlignmentReader::ReaderCount << std::endl;
				// std::cout << contigLength << std::endl;
				// ThreadPool::Instance()->joinAll();
				// std::cout << "starting back" << std::endl;
				// ThreadPool::Instance()->startIOService();
			// }
			if (count % 100000 == 0)
			{
				std::cout << count << std::endl;
			}
			// if there is a gap between the variants that is larger than the padding then generate a graph and push it in the queue
			if ((previousVariantEndPosition > 0) && (contigLength > minContig) && (this->m_padding < variantPositionDifference))
			{
				/*
				if (count == 86458)
				{
					std::cout << "inside" << std::endl;
					std::cout << contigLength << std::endl;
				}
				*/
				// std::cout << count << " " << contigLength << std::endl;

				auto alignmentReaderPtr = this->m_alignment_reader_manager->generateAlignmentReader(); // create alignment reader

				// create region for alignmentReader
				std::string startPosition = std::to_string(firstVariantInList->getPosition() - this->m_padding);
				std::string endPosition = std::to_string(previousVariantEndPosition);
				auto region = std::make_shared< Region >(this->m_reference_ptr->getRegion()->getReferenceID() + ":" + startPosition + "-" + endPosition);

				alignmentReaderPtr->setRegion(region); // set alignmentreader's region
				auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantListPtr, alignmentReaderPtr);
				/*
				if (count >= 98957)
				{
					std::cout << "WOW: " << count << std::endl;
					ThreadPool::Instance()->m_print_stuff = true;
				}
				*/
				// this->m_gssw_graphs.push(gsswGraph);
				ThreadPool::Instance()->postJob(boost::bind(&GSSWGraph::constructGraph, gsswGraph));
				// ThreadPool::Instance()->postJob(boost::bind(&IAlignmentReader::releaseReader, alignmentReaderPtr));
				variantListPtr = std::make_shared< VariantList >(); // contains the variants for contiguous sections of variants
				firstVariantInList = variantPtr; // set the first variant
				contigLength = 0;
			}
			variantListPtr->addVariant(variantPtr);
			previousVariantEndPosition = variantPtr->getPosition() + variantPtr->getRef()[0].size();
		}
		std::cout << "finished" << std::endl;
		ThreadPool::Instance()->joinAll();
		std::cout << "threads finished" << std::endl;
	}
}
}
