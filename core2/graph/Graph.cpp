#include "Graph.h"

#include "core2/util/Types.h"

namespace graphite
{
	Graph::Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing) :
		m_fasta_reference_ptr(fastaReferencePtr),
		m_variant_ptrs(variantPtrs),
		m_graph_spacing(graphSpacing)
	{
	}

	Graph::~Graph()
	{
	}

	void Graph::generateGraph()
	{
		std::string referenceSequence;
		Region::SharedPtr referenceRegionPtr;
		getGraphReference(referenceSequence, referenceRegionPtr);
		Node::SharedPtr firstNodePtr;
		Node::SharedPtr lastNodePtr;
		generateReferenceGraphNode(firstNodePtr, lastNodePtr, referenceSequence, referenceRegionPtr);
		addVariantsToGraph(firstNodePtr);
		firstNode = condenseGraph(lastNodePtr);
		// calculate prefix and suffix matching sequences
		// add gssw mechanism
	}

	std::vector< Region::SharedPtr > Graph::getRegionPtrs()
	{
		return this->m_graph_regions;
	}

	void Graph::adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr)
	{
	}

	void Graph::getGraphReference(std::string& sequence, Region::SharedPtr regionPtr)
	{
		std::string referenceID = this->m_variant_ptrs[0]->getChromosome();
		position startPosition = MAX_POSITION;
		position endPosition = 0;
		for (auto variantPtr : this->m_variant_ptrs)
		{
			if (variantPtr->getPosition() < startPosition)
			{
				startPosition = variantPtr->getPosition();
			}
			position variantEndPosition = variantPtr->getPosition() + variantPtr->getReferenceAllelePtr()->getSequence().size();
			if (endPosition < variantEndPosition)
			{
				endPosition = variantEndPosition;
			}
		}
		startPosition -= this->m_graph_spacing;
		endPosition += this->m_graph_spacing;
		regionPtr = std::make_shared< Region >(referenceID, startPosition, endPosition, Region::BASED::ONE);
		sequence = this->m_fasta_reference_ptr->getSequenceStringFromRegion(regionPtr);
	}

	void Graph::generateReferenceGraphNode(Node::SharedPtr firstNodePtr, Node::SharedPtr lastNodePtr, const std::string& referenceSequence, Region::SharedPtr regionPtr)
	{
		firstNodePtr = nullptr;
		position nodePosition = regionPtr->getStartPosition();
		Node::SharedPtr prevNode = nullptr;
		for (auto i = 0; i < referenceSequence.size(); ++i)
		{
			Node::SharedPtr node = std::make_shared< Node >(referenceSequence.c_str()[i], 1, nodePosition, Node::ALLELE_TYPE::REF);
			if (firstNodePtr == nullptr)
			{
				firstNodePtr = node;
			}
			if (prevNode != nullptr)
			{
				node->addInNode(prevNode);
				prevNode->addOutNode(node);
			}
			prevNode = node;
			nodePosition += 1;
			lastNodePtr = node;
		}
	}

	void Graph::addVariantsToGraph(Node::SharedPtr firstNodePtr)
	{
		std::unordered_map< position, Node::SharedPtr > referenceNodePtrPositionMap;
		Node::SharedPtr nodePtr = firstNodePtr;
		do
		{
			referenceNodePtrPositionMap.emplace(nodePtr->getPosition(), nodePtr);
			if (nodePtr->getOutNodes().size() > 0)
			{
				nodePtr = nodePtr->getOutNodes()[0];
			}
		} while (nodePtr->getOutNodes().size() > 0);


		for (auto variantPtr : this->m_variant_ptrs)
		{
			position variantPosition = variantPtr->getPosition();
			position inPosition = variantPtr->getReferenceAllelePtr()->getSequence().size() + 1;
			for (uint32_t i = variantPosition; i < inPosition; ++i)
			{
				auto iter = referenceNodePtrPositionMap.find(i);
				if (iter != referenceNodePtrPositionMap.end())
				{
					iter->second->addOverlappingAllelePtr(variantPtr->getReferenceAllelePtr());
				}
			}
			auto inReferenceIter = referenceNodePtrPositionMap.find(variantPosition);
			auto outReferenceIter = referenceNodePtrPositionMap.find(inPosition);
			if (inReferenceIter == referenceNodePtrPositionMap.end() || outReferenceIter == referenceNodePtrPositionMap.end())
			{
				std::cout << "Invalid Graph: addVariantsToGraph, position: " << variantPtr->getPosition() << std::endl;
				exit(EXIT_FAILURE);
			}
			for (auto altAllelePtr : variantPtr->getAlternateAllelePtrs())
			{
				auto altNodePtr = std::make_shared< Node >(altAllelePtr->getSequence(), variantPosition, Node::ALLELE_TYPE::ALT);
				altNodePtr->addOverlappingAllelePtr(variantPtr->getReferenceAllelePtr());
				(inReferenceIter->second)->addOutNode(altNodePtr);
				(outReferenceIter->second)->addInNode(altNodePtr);
			}
		}
	}
	Node::SharedPtr Graph::condenseGraph(Node::SharedPtr lastNodePtr)
	{
		Node::SharedPtr nodePtr = lastNodePtr;
		do
		{
			Node::SharedPtr refNodePtr = nodePtr->getReferenceInNode();
			if (refNodePtr == nullptr)
			{
				std::cout << "Invalid Graph: condenseGraph, position: " << nodePtr->getPosition() << std::endl;
				exit(EXIT_FAILURE);
			}
			nodePtr = Node::mergeNodes(refNodePtr, nodePtr);
		} while (nodePtr->getInNodes().size() > 0);
		return nodePtr;
	}
}
