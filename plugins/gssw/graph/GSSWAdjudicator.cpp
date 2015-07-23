#include "GSSWAdjudicator.h"
#include "GSSWGraph.h"
#include "core/variant/VariantList.h"
#include "AlignmentReporter.h"

#include <memory>

namespace gwiz
{
namespace gssw
{
	GSSWAdjudicator::GSSWAdjudicator(uint32_t swPercent) :
		m_sw_percent(swPercent)
	{
		static int id = 0;
		m_id = ++id;
	}

	GSSWAdjudicator::~GSSWAdjudicator()
	{
	}

	void GSSWAdjudicator::printNodes(GSSWGraph::SharedPtr graphPtr, const std::string& alignment)
	{
		auto gsswGraph = graphPtr->getGSSWGraph();
		for (uint32_t i = 0; i < gsswGraph->size; ++i)
		{
			auto node = gsswGraph->nodes[i];
			std::cout << "pos: " << node->position << std::endl;
			gssw_print_score_matrix(node->seq, node->len, alignment.c_str(), alignment.size(), node->alignment);
		}
	}

	IVariantList::SharedPtr GSSWAdjudicator::adjudicateGraph(IGraph::SharedPtr graphPtr, IAlignmentList::SharedPtr alignmentListPtr)
	{
		throw "GSSWAdjudicator::adjudicateGraph needs to be fixed";
		auto gsswGraphPtr = std::dynamic_pointer_cast< GSSWGraph >(graphPtr);
		if (gsswGraphPtr) // kind of punting for now. in the future this should be updated so it handles all igraphs the same
		{
			IAlignment::SharedPtr alignmentPtr;
			while (alignmentListPtr->getNextAlignment(alignmentPtr))
			{
				auto graphMappingPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr);
				gssw_node_cigar* nc = graphMappingPtr->cigar.elements;
				// printNodes(gsswGraphPtr, std::string(alignmentPtr->getSequence(), alignmentPtr->getLength()));
				float percentage = (this->m_sw_percent * 0.01);
				bool mapped = (graphMappingPtr->score >= ((alignmentPtr->getLength() * gsswGraphPtr->getMatchValue()) * percentage));
				std::vector< std::tuple< uint32_t, std::string > > variantInformation;
				for (int i = 0; i < graphMappingPtr->cigar.length; ++i, ++nc)
				{
					IAllele::SharedPtr allelePtr = *(std::shared_ptr< IAllele >*)(nc->node->data);
					std::cout << allelePtr->getSequence() << std::endl;
					/*
					auto variantPtr = gsswGraphPtr->getVariantFromNodeID(nc->node->id);
					if (variantPtr != nullptr)
					{
						if (mapped)
						{
							variantInformation.emplace_back(std::make_tuple< uint32_t, std::string >(variantPtr->getVariantID(), std::string(nc->node->seq, nc->node->len)));
						}
						variantPtr->addPotentialAlignment(alignmentPtr);
					}
					*/
				}
				if (mapped)
				{
					alignmentPtr->setMappingInformation(graphMappingPtr->score, variantInformation);
				}
			}
		}
		else
		{
			throw "adjudicateGraph has not been implemented for non-GSSWGraphs";
		}
		return nullptr;
	}
}
}
