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
		// std::mutex reportedVariantsMutex;

		std::string referenceID = regionPtr->getReferenceID();
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();
		position currentPosition = startPosition;

		while (currentPosition < endPosition)
		{
			auto endGraphPosition = (currentPosition + graphSize > endPosition) ? endPosition : (currentPosition + graphSize);
			auto graphRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition) + "-" + std::to_string(endGraphPosition)));
			auto variantsListPtr = this->m_variant_list_ptr->getVariantsInRegion(graphRegion);
			if (variantsListPtr->getCount() > 0) // if we have variants, then process them
			{
				/*
				static uint32_t counter = 0;
				if (++counter % 1000 == 0)
				{
					std::cout << "count: " << counter << std::endl;
				}
				*/
				// std::cout << "asdf1" << std::endl;

				auto alignmentReaderPtr = this->m_alignment_reader_manager->generateAlignmentReader(); // create alignment reader
				auto alignmentRegion = std::make_shared< Region >(std::string(referenceID + ":" + std::to_string(currentPosition + alignmentPadding) + "-" + std::to_string(endGraphPosition - alignmentPadding)));
				// std::cout << "asdf2" << std::endl;

				// create region for alignmentReader
				alignmentReaderPtr->init();
				// std::cout << "asdf3" << std::endl;
				alignmentReaderPtr->setRegion(alignmentRegion); // set alignmentreader's region
				if (alignmentReaderPtr->getReadCount() == 0) { continue; } // if there are no reads then continue

				// std::cout << "asdf4" << std::endl;

				auto funct = std::bind(&GraphManager::constructAndAdjudicateGraph, this, reportedVariants, variantsListPtr, alignmentReaderPtr, currentPosition, graphSize);
				ThreadPool::Instance()->enqueue(funct);

				// std::cout << "asdf5" << std::endl;

			}
			currentPosition += graphSize - overlap;
		}
		ThreadPool::Instance()->joinAll();
		reportedVariants->sort();
		return reportedVariants;
	}

	void GraphManager::constructAndAdjudicateGraph(IVariantList::SharedPtr reportedVariants, IVariantList::SharedPtr variantsListPtr, IAlignmentReader::SharedPtr alignmentReaderPtr, position startPosition, size_t graphSize)
	{
		auto gsswGraph = std::make_shared< GSSWGraph >(this->m_reference_ptr, variantsListPtr, startPosition, graphSize);
		gsswGraph->constructGraph();

		auto variantsPtrList =  this->m_graph_adjudicator_ptr->adjudicateGraph(gsswGraph, alignmentReaderPtr);

		{
			this->m_gssw_graph_mutex.lock();
			reportedVariants->addVariants(variantsPtrList);
			this->m_gssw_graph_mutex.unlock();
		}
	}

}
}
