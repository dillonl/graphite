#include "Graph.h"

#include "core2/util/Types.h"

namespace graphite
{
	Graph::Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing) :
		m_fasta_reference_ptr(fastaReferencePtr),
		m_variant_ptrs(variantPtrs),
		m_graph_spacing(graphSpacing),
		m_score_threshold(70)
	{
		generateGraph();
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
		firstNodePtr = condenseGraph(lastNodePtr);
		setPrefixAndSuffix(firstNodePtr); // calculate prefix and suffix matching sequences
		this->m_first_node = firstNodePtr;
	}

	std::vector< Region::SharedPtr > Graph::getRegionPtrs()
	{
		return this->m_graph_regions;
	}

	void Graph::getGraphReference(std::string& sequence, Region::SharedPtr& regionPtr)
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

	void Graph::generateReferenceGraphNode(Node::SharedPtr& firstNodePtr, Node::SharedPtr& lastNodePtr, const std::string& referenceSequence, Region::SharedPtr regionPtr)
	{
		firstNodePtr = nullptr;
		position nodePosition = regionPtr->getStartPosition();
		Node::SharedPtr prevNode = nullptr;
		for (auto i = 0; i < referenceSequence.size(); ++i)
		{
			Node::SharedPtr node = std::make_shared< Node >(referenceSequence.c_str() + i, 1, nodePosition, Node::ALLELE_TYPE::REF);
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
				nodePtr = *nodePtr->getOutNodes().begin();
			}
		} while (nodePtr->getOutNodes().size() > 0);


		for (auto variantPtr : this->m_variant_ptrs)
		{
			position variantPosition = variantPtr->getPosition();
			position inPosition = variantPosition + variantPtr->getReferenceAllelePtr()->getSequence().size() + 1;
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
				altAllelePtr->registerNodePtr(altNodePtr);
				altNodePtr->addOverlappingAllelePtr(variantPtr->getReferenceAllelePtr());
				altNodePtr->addInNode(inReferenceIter->second);
				altNodePtr->addOutNode(outReferenceIter->second);
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
			if ((nodePtr->getInNodes().size() > 1) || (refNodePtr->getOutNodes().size() > 1))
			{
				nodePtr = refNodePtr;
			}
			else
			{
				nodePtr = Node::mergeNodes(refNodePtr, nodePtr);
			}
		} while (nodePtr->getInNodes().size() > 0);
		return nodePtr;
	}

	void Graph::setPrefixAndSuffix(Node::SharedPtr firstNodePtr)
	{
		Node::SharedPtr nextRefPtr = firstNodePtr;
		this->m_node_ptrs_map.emplace(firstNodePtr->getID(), firstNodePtr);
		std::unordered_set< Node::SharedPtr > outNodePtrs = nextRefPtr->getOutNodes();
		while (outNodePtrs.size() > 0)
		{
			for (auto node1Ptr : outNodePtrs)
			{
				this->m_node_ptrs_map.emplace(node1Ptr->getID(), node1Ptr);
				if (node1Ptr->getAlleleType() == Node::ALLELE_TYPE::REF)
				{
					nextRefPtr = node1Ptr;
					if (outNodePtrs.size() == 1) // skip comparion to self
					{
						continue;
					}
				}
				for (auto node2Ptr : outNodePtrs)
				{
					if (node1Ptr == node2Ptr)
					{
						continue;
					}
					const char* seq1;
					const char* seq2;
					size_t len1;
					size_t len2;
					if (node1Ptr->getSequence().size() < node2Ptr->getSequence().size())
					{
						seq1 = node1Ptr->getSequence().c_str();
						seq2 = node2Ptr->getSequence().c_str();
						len1 = node1Ptr->getSequence().size();
						len2 = node2Ptr->getSequence().size();
					}
					else
					{
						seq2 = node1Ptr->getSequence().c_str();
						seq1 = node2Ptr->getSequence().c_str();
						len2 = node1Ptr->getSequence().size();
						len1 = node2Ptr->getSequence().size();
					}
					// seq1 and len1 is less than seq2 and len2
					size_t maxPrefix;
					for (maxPrefix = 0; maxPrefix < len1; ++maxPrefix)
					{
						if (seq1[maxPrefix] != seq2[maxPrefix])
						{
							break;
						}
					}
					if (node1Ptr->getIdenticalPrefixLength() < maxPrefix)
					{
						node1Ptr->setIdenticalPrefixLength(maxPrefix);
					}
					if (node2Ptr->getIdenticalPrefixLength() < maxPrefix)
					{
						node2Ptr->setIdenticalPrefixLength(maxPrefix);
					}
					size_t maxSuffix;
					for (maxSuffix = (len1 - 1); maxSuffix >= 0; --maxSuffix)
					{
						if (seq1[maxSuffix] != seq2[maxSuffix])
						{
							break;
						}
					}
					if (node1Ptr->getIdenticalSuffixLength() < maxSuffix)
					{
						node1Ptr->setIdenticalSuffixLength(maxSuffix);
					}
					if (node2Ptr->getIdenticalSuffixLength() < maxSuffix)
					{
						node2Ptr->setIdenticalSuffixLength(maxSuffix);
					}
				}
			}
			outNodePtrs = nextRefPtr->getOutNodes();
		}
	}

	void Graph::adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue)
	{
		int8_t* nt_table = gssw_create_nt_table();
		int8_t* mat = gssw_create_score_matrix(matchValue, mismatchValue);
		gssw_graph* graph = gssw_graph_create(this->m_node_ptrs_map.size());

		unordered_map< uint32_t, gssw_node* > gsswNodePtrsMap;
		for (auto iter : this->m_node_ptrs_map)
		{
			Node::SharedPtr nodePtr = iter.second;
			gssw_node* gsswNode = (gssw_node*)gssw_node_create(nodePtr.get(), 0, nodePtr->getSequence().c_str(), nullptr, nullptr);
			gsswNodePtrsMap.emplace(nodePtr->getID(), gsswNode);
			gssw_graph_add_node(graph, gsswNode);
		}

		for (auto iter : this->m_node_ptrs_map)
		{
			Node::SharedPtr nodePtr = iter.second;
			gssw_node* gsswNode = gsswNodePtrsMap[nodePtr->getID()];
			for (auto inNodePtr : nodePtr->getInNodes())
			{
				gssw_node* gsswInNode = gsswNodePtrsMap[inNodePtr->getID()];
				gssw_nodes_add_edge(gsswInNode, gsswNode);
			}
			for (auto outNodePtr : nodePtr->getOutNodes())
			{
				gssw_node* gsswOutNode = gsswNodePtrsMap[outNodePtr->getID()];
				gssw_nodes_add_edge(gsswNode, gsswOutNode);
			}
		}

		gssw_graph_fill(graph, bamAlignmentPtr->QueryBases.c_str(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0, 15, 2, true);
		gssw_graph_mapping* gm = gssw_graph_trace_back (graph, bamAlignmentPtr->QueryBases.c_str(), bamAlignmentPtr->QueryBases.size(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0);
		processTraceback(gm, bamAlignmentPtr, !bamAlignmentPtr->IsReverseStrand(), matchValue, mismatchValue, gapOpenValue, gapExtensionValue);

		gssw_graph_mapping_destroy(gm);

		// note that nodes which are referred to in this graph are destroyed as well
		gssw_graph_destroy(graph);

		free(nt_table);
		free(mat);
	}

	void Graph::processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue)
	{
		std::vector< std::tuple< Node*, uint32_t > > nodePtrScoreTuples;
		uint32_t prefixMatch = 0;
		uint32_t suffixMatch = 0;
		bool setPrefix = true;
		uint32_t totalScore = 0;
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			gssw_node* gsswNode = nc[i].node;
			Node* nodePtr = (Node*)gsswNode->data;
			int32_t score = 0;
			uint32_t length = 0;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					suffixMatch += nc->cigar->elements[j].length;
					if (setPrefix)
					{
						prefixMatch += nc->cigar->elements[j].length;
					}
					score += (matchValue * nc->cigar->elements[j].length);
					break;
				case 'S':
				case 'X':
					setPrefix = false;
					suffixMatch = 0;
					score -= (mismatchValue * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					setPrefix = false;
					suffixMatch = 0;
					score -= gapOpenValue;
					score -= (gapExtensionValue * (nc->cigar->elements[j].length -1));
					break;
				default:
					break;
				}
				length += nc->cigar->elements[j].length;
			}
			score = (score < 0) ? 0 : score; // the floor of the mapping score is 0
			std::tuple< Node*, uint32_t > nodePtrScoreTuple(nodePtr, score);
			nodePtrScoreTuples.emplace_back(nodePtrScoreTuple);
			totalScore += score;
		}

		for (auto nodePtrScoreTuple : nodePtrScoreTuples)
		{
			Node* nodePtr = std::get< 0 >(nodePtrScoreTuple);
			uint32_t nodeScore = std::get< 1 >(nodePtrScoreTuple);
			if (nodePtr->getIdenticalPrefixLength() <= prefixMatch || nodePtr->getIdenticalSuffixLength() <= suffixMatch)
			{
				nodePtr->incrementScoreCount(bamAlignmentPtr, isForwardStrand, -1);
			}
			else if (totalScore < m_score_threshold)
			{
				nodePtr->incrementScoreCount(bamAlignmentPtr, isForwardStrand, 0);
			}
			else
			{
				nodePtr->incrementScoreCount(bamAlignmentPtr, isForwardStrand, nodeScore);
			}
		}
	}

}
