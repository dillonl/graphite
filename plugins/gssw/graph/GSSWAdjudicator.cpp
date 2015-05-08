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

	IVariantList::SharedPtr GSSWAdjudicator::adjudicateGraph(IGraph::SharedPtr graphPtr, IAlignmentReader::SharedPtr alignmentsReaderPtr)
	{
		// ADD AN ID TO TRACK WHEN THIS FALLS OUT OF SCOPE
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
				bool mapped = (graphMappingPtr->score >= ((alignmentPtr->getLength() * gsswGraphPtr->getMatchValue()) * 0.75));
				std::unordered_map< uint32_t, int32_t > alignmentIDVariants;
				for (int i = 0; i < graphMappingPtr->cigar.length; ++i, ++nc)
				{
					auto variantPtr = gsswGraphPtr->getVariantFromNodeID(nc->node->id);
					if (variantPtr != nullptr)
					{
						this->m_adjudication_lock.lock();
						alignmentPtr->setMapped(mapped);
						alignmentIDVariants.emplace(std::make_pair(variantPtr->getVariantID(), graphMappingPtr->score));
						variantPtr->addPotentialAlignment(alignmentPtr, std::string(nc->node->seq, nc->node->len));
						variantList->addVariant(variantPtr);
						this->m_adjudication_lock.unlock();
					}
				}
				{
					// std::cout << "locking: " << m_id << std::endl;
					this->m_adjudication_lock.lock();
					alignmentPtr->setMappingScoreAndVariants(graphMappingPtr->score, alignmentIDVariants);
					this->m_adjudication_lock.unlock();
					// std::cout << "unlocking: " << m_id << std::endl;
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
