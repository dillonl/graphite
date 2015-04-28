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
		auto variantList = std::make_shared< VariantList >();
		auto gsswGraphPtr = std::dynamic_pointer_cast< GSSWGraph >(graphPtr);
		if (gsswGraphPtr) // kind of punting for now. in the future this should be updated so it handles all igraphs the same
		{
			IAlignment::SharedPtr alignmentPtr;
			while (alignmentsReaderPtr->getNextAlignment(alignmentPtr))
			{
				auto graphMappingPtr = gsswGraphPtr->traceBackAlignment(alignmentPtr);
				gssw_node_cigar* nc = graphMappingPtr->cigar.elements;
				// printNodes(gsswGraphPtr, std::string(alignmentPtr->getSequence(), alignmentPtr->getLength()));
				bool pass = (graphMappingPtr->score >= ((alignmentPtr->getLength() * gsswGraphPtr->getMatchValue()) * 0.75));
				for (int i = 0; i < graphMappingPtr->cigar.length; ++i, ++nc)
				{
					auto variantPtr = gsswGraphPtr->getVariantFromNodeID(nc->node->id);
					if (variantPtr != nullptr)
					{
						if (pass && !variantPtr->getPass()) // since this variant is adjudicated by multiple alignments set it to pass only if the variants pass is false (we do not want to set a passing variant to false)
						{
							variantPtr->setPass(true);
						}
						this->m_adjudication_lock.lock();
						if (pass) // only count the variants tha pass
						{
							variantPtr->increaseCount(nc->node->seq, nc->node->len, alignmentPtr); // record the variant (ref or alt) that went through the node
						}
						variantList->addVariant(variantPtr);
						this->m_adjudication_lock.unlock();
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
