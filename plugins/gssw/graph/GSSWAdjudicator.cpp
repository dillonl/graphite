#include "GSSWAdjudicator.h"
#include "GSSWGraph.h"
#include "core/variants/VariantList.h"
#include "AlignmentReporter.h"

namespace gwiz
{
namespace gssw
{
	GSSWAdjudicator::GSSWAdjudicator()
	{
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

	IVariantList::SharedPtr GSSWAdjudicator::adjudicateGraph(IGraph::SharedPtr graphPtr, IAlignmentReader::SharedPtr alignmentsReaderPtr)
	{
		static std::mutex adjLock;

		auto variantList = std::make_shared< VariantList >();
		auto gsswGraphPtr = std::dynamic_pointer_cast< GSSWGraph >(graphPtr);
		// std::cout << "adj" << std::endl;
		if (gsswGraphPtr) // kind of punting for now. in the future this should be updated so it handles all igraphs the same
		{
			// std::cout << "adj2" << std::endl;
			IAlignment::SharedPtr alignmentPtr;
			while (alignmentsReaderPtr->getNextAlignment(alignmentPtr))
			{
				auto graphMappingPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr);
				// gsswGraphPtr->recordAlignmentVariants(graphMappingPtr, alignmentPtr);
				gssw_node_cigar* nc = graphMappingPtr->cigar.elements;
				// printNodes(gsswGraphPtr, std::string(alignmentPtr->getSequence(), alignmentPtr->getLength()));
				if (graphMappingPtr->score < ((alignmentPtr->getLength() * gsswGraphPtr->getMatchValue()) * 0.75)) // skip low scoring	mappings
				{
					continue;
				}
				// std::cout << graphMappingPtr->score << " " <<   (this->m_max_mapping_score * 0.9) << std::endl;
				for (int i = 0; i < graphMappingPtr->cigar.length; ++i, ++nc)
				{
					auto variantPtr = gsswGraphPtr->getVariantFromNodeID(nc->node->id);
					if (variantPtr != nullptr)
					{
						std::unique_lock< std::mutex > lock(adjLock);
						variantPtr->increaseCount(nc->node->seq, nc->node->len, alignmentPtr); // record the variant (ref or alt) that went through the node
						variantList->addVariant(variantPtr);
					}
				}
			}
		}
		else
		{
			throw "adjudicateGraph has not been implemented for non-GSSWGraphs";
		}
		// AlignmentReporter::Instance()->printAlignmentReportsToStream(std::cout);
		return variantList;
	}
}
}
