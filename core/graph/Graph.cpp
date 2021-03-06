#include "Graph.h"
// #include "Traceback.h"
#include "GraphTraceback.hpp"

#include "core/util/Types.h"
#include <deque>

namespace graphite
{
	Graph::Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing, bool printGraph) :
		m_fasta_reference_ptr(fastaReferencePtr),
		// m_variant_ptrs(variantPtrs),
		m_graph_spacing(graphSpacing),
		m_score_threshold(70),
		m_graph_printer_ptr(nullptr)
	{
		m_variant_ptrs = reconcileVariantSemantics(variantPtrs);
		generateGraph();
		if (printGraph)
		{
			m_graph_printer_ptr = std::make_shared< GraphPrinter >(this);
		}
	}

	Graph::~Graph()
	{
		if (m_graph_printer_ptr != nullptr)
		{
			m_graph_printer_ptr->printGraph();
		}
		m_fasta_reference_ptr = nullptr;
		m_first_node = nullptr;
		m_variant_ptrs.clear();
		m_graph_regions.clear();
		m_node_ptrs_map.clear();
		for (auto nodePtr : m_all_created_nodes)
		{
			nodePtr->clearInAndOutNodes();
		}
		m_all_created_nodes.clear();
	}

	void Graph::generateGraph()
	{
		std::string referenceSequence;
		Region::SharedPtr referenceRegionPtr;
		getGraphReference(referenceSequence, referenceRegionPtr, this->m_variant_ptrs);
		Node::SharedPtr firstNodePtr;
		Node::SharedPtr lastNodePtr;
		generateReferenceGraphNode(firstNodePtr, lastNodePtr, referenceSequence, referenceRegionPtr);
		addVariantsToGraph(firstNodePtr);
		firstNodePtr = condenseGraph(lastNodePtr);
		setPrefixAndSuffix(firstNodePtr); // calculate prefix and suffix matching sequences
		this->m_first_node = firstNodePtr;
		setRegionPtrs();
		addAdditionalNodeEdges();
	}

	std::vector< Variant::SharedPtr > Graph::reconcileVariantSemantics(std::vector< Variant::SharedPtr >& variantPtrs)
	{
		std::vector< Variant::SharedPtr > uniqueVariantPtrs;
		std::string referenceSequence;
		Region::SharedPtr referenceRegionPtr;
		getGraphReference(referenceSequence, referenceRegionPtr, variantPtrs);

		std::unordered_map< std::string, Allele::SharedPtr > variantSequenceMap;
		std::unordered_map< std::string, Variant::SharedPtr > variantPtrMap;
		std::unordered_set< Variant::SharedPtr > variantPtrSet;
		// variantSequences.emplace_back(referenceSequence); // we aren't putting in the reference allele because there are many ref alleles and they all have the same equence so they will be over counted
		for (auto variantPtr : variantPtrs)
		{
			int prefixSize = variantPtr->getPosition() - referenceRegionPtr->getStartPosition();
			int suffixSize = (prefixSize + variantPtr->getReferenceAllelePtr()->getSequence().size());
			std::string prefix = referenceSequence.substr(0, prefixSize);
			std::string suffix = referenceSequence.substr(suffixSize);
			Allele::SharedPtr refAllelePtr = variantPtr->getReferenceAllelePtr();
			for (auto altAllelePtr : variantPtr->getAlternateAllelePtrs())
			{
				std::string sequence = prefix + altAllelePtr->getSequence() + suffix;
				auto variantIter = variantSequenceMap.find(sequence);
				if (variantIter != variantSequenceMap.end()) // if there is a dup sequence that means there is a semantic sequence match so pair the alleles
				{
					auto variantPtrIter = variantPtrMap.find(sequence);
					if (variantPtrIter == variantPtrMap.end()) { std::cout << "there was a problem in Graph.cpp." << std::endl; }
					Variant::SharedPtr queryVariantPtr = variantPtrIter->second;
					Allele::SharedPtr allelePtr = variantIter->second;
					Allele::SharedPtr queryRefAllelePtr = queryVariantPtr->getReferenceAllelePtr();
					queryRefAllelePtr->pairAllele(refAllelePtr); // pair the ref allele so they are both counted when there is a SW match
					allelePtr->pairAllele(altAllelePtr); // pair the allele so they are both counted when there is a SW match
					allelePtr->addSemanticLoci(variantPtr->getPosition(), refAllelePtr->getSequence(), altAllelePtr->getSequence());
					altAllelePtr->addSemanticLoci(queryVariantPtr->getPosition(), queryRefAllelePtr->getSequence(), allelePtr->getSequence());
				}
				else
				{
					variantPtrMap.emplace(sequence, variantPtr);
					variantSequenceMap.emplace(sequence, altAllelePtr);
					if (variantPtrSet.find(variantPtr) == variantPtrSet.end()) // we don't want to add the variant more than once on multi-allelics
					{
						uniqueVariantPtrs.emplace_back(variantPtr);
					}
				}
			}
		}
		return uniqueVariantPtrs;
	}

	std::vector< Region::SharedPtr > Graph::getRegionPtrs()
	{
		return this->m_graph_regions;
	}

	void Graph::compressLargeNodes()
	{
		std::deque< Node::SharedPtr > nodes;
		nodes.emplace_back(this->m_first_node);
		while (!nodes.empty())
		{
			Node::SharedPtr nodePtr = nodes.front();
			nodes.pop_front();
			if (nodePtr->getSequence().size() > (this->m_graph_spacing * 2))
			{
				std::string newSequence = nodePtr->getSequence().substr(0, this->m_graph_spacing);
				newSequence += std::string(this->m_graph_spacing, 'N');
				newSequence += nodePtr->getSequence().substr(this->m_graph_spacing);
				nodePtr->setCompressedSequence(newSequence);
			}
			for (auto outNodePtr : nodePtr->getOutNodes())
			{
				nodes.emplace_back(outNodePtr);
			}
		}
	}

	void Graph::setRegionPtrs()
	{
		this->m_graph_regions.clear();
		std::string referenceID = this->m_variant_ptrs[0]->getChromosome();
		Node::SharedPtr currentRefNode = this->m_first_node;
		position startPosition = currentRefNode->getPosition();
		position endPosition = currentRefNode->getPosition() + currentRefNode->getOriginalSequenceSize();
		Region::SharedPtr currentRegionPtr = nullptr;
		while (currentRefNode != nullptr)
		{
			if (currentRefNode->getSequence().size() != currentRefNode->getOriginalSequenceSize())
			{
				endPosition = currentRefNode->getPosition() + this->m_graph_spacing;
				this->m_graph_regions.emplace_back(std::make_shared< Region >(referenceID, startPosition, endPosition, Region::BASED::ONE));
				startPosition = currentRefNode->getPosition() +  currentRefNode->getOriginalSequenceSize() - this->m_graph_spacing;
			}
			else
			{
				endPosition = currentRefNode->getPosition() + currentRefNode->getSequence().size();
			}
			currentRefNode = currentRefNode->getReferenceOutNode();
		}
		this->m_graph_regions.emplace_back(std::make_shared< Region >(referenceID, startPosition, endPosition, Region::BASED::ONE));
	}

	void Graph::getGraphReference(std::string& sequence, Region::SharedPtr& regionPtr, std::vector< Variant::SharedPtr >& variantPtrs)
	{
		std::string referenceID = variantPtrs[0]->getChromosome();
		position startPosition = MAX_POSITION;
		position endPosition = 0;
		for (auto variantPtr : variantPtrs)
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
		int tmpStartPos = startPosition - this->m_graph_spacing;
		startPosition = (tmpStartPos < 0) ? 1 : startPosition -= this->m_graph_spacing;
		endPosition += this->m_graph_spacing;
		regionPtr = std::make_shared< Region >(referenceID, startPosition, endPosition, Region::BASED::ONE);
		sequence = this->m_fasta_reference_ptr->getSequenceStringFromRegion(regionPtr);
	}

	void Graph::generateReferenceGraphNode(Node::SharedPtr& firstNodePtr, Node::SharedPtr& lastNodePtr, const std::string& referenceSequence, Region::SharedPtr regionPtr)
	{
		firstNodePtr = nullptr;
		position nodePosition = regionPtr->getStartPosition();
		Node::SharedPtr prevNodePtr = nullptr;
		static bool firstTime = false;
		for (auto i = 0; i < referenceSequence.size(); ++i)
		{
			Node::SharedPtr nodePtr = std::make_shared< Node >(referenceSequence.c_str() + i, 1, nodePosition, Node::ALLELE_TYPE::REF);
			this->m_all_created_nodes.emplace(nodePtr);
			if (firstNodePtr == nullptr)
			{
				firstNodePtr = nodePtr;
			}
			if (prevNodePtr != nullptr)
			{
				nodePtr->addInNode(prevNodePtr);
				prevNodePtr->addOutNode(nodePtr);
			}
			prevNodePtr = nodePtr;
			nodePosition += 1;
			lastNodePtr = nodePtr;
		}
		firstTime = false;
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

		// std::vector< Node::SharedPtr > previousNodes;
		std::unordered_map< position, std::vector< Node::SharedPtr > > endPositionNodePtrs;
		for (auto variantPtr : this->m_variant_ptrs)
		{
			position variantNodeStartPosition = variantPtr->getPosition();
			position variantNodeEndPosition = variantPtr->getPosition() + variantPtr->getReferenceAllelePtr()->getSequence().size() - 1;
			position leftAdjacentNodeEndPosition = variantNodeStartPosition - 1;
			position rightAdjacentNodeStartPosition = variantNodeEndPosition + 1;
			for (uint32_t i = variantNodeStartPosition; i <= variantNodeEndPosition; ++i)
			{
				auto iter = referenceNodePtrPositionMap.find(i);
				if (iter != referenceNodePtrPositionMap.end())
				{
					iter->second->registerAllelePtr(variantPtr->getReferenceAllelePtr());
				}
			}
			auto leftAdjacentReferenceNodeIter = referenceNodePtrPositionMap.find(leftAdjacentNodeEndPosition);
			auto rightAdjacentReferenceNodeIter = referenceNodePtrPositionMap.find(rightAdjacentNodeStartPosition);
			if (leftAdjacentReferenceNodeIter == referenceNodePtrPositionMap.end() || rightAdjacentReferenceNodeIter == referenceNodePtrPositionMap.end())
			{
				std::cout << "Invalid Graph: addVariantsToGraph, position: " << leftAdjacentNodeEndPosition << " - "  << rightAdjacentNodeStartPosition << std::endl;
				exit(EXIT_FAILURE);
			}
			for (auto altAllelePtr : variantPtr->getAlternateAllelePtrs())
			{
				bool duplicateVariantFound = false;
				for (auto siblingNodePtr : leftAdjacentReferenceNodeIter->second->getOutNodes())
				{
					if (siblingNodePtr->getAlleleType() == Node::ALLELE_TYPE::ALT && siblingNodePtr->getSequence().compare(altAllelePtr->getSequence()) == 0) // if we find an exact match then just register the node and move on
					{
						duplicateVariantFound = true;
						siblingNodePtr->registerAllelePtr(altAllelePtr);
						altAllelePtr->registerNodePtr(siblingNodePtr);
						break;
					}
				}
				if (!duplicateVariantFound)
				{
					auto altNodePtr = std::make_shared< Node >(altAllelePtr->getSequence(), variantNodeStartPosition, Node::ALLELE_TYPE::ALT);
					this->m_all_created_nodes.emplace(altNodePtr);
					altNodePtr->registerAllelePtr(altAllelePtr);
					altAllelePtr->registerNodePtr(altNodePtr);

					altNodePtr->addInNode(leftAdjacentReferenceNodeIter->second);
					altNodePtr->addOutNode(rightAdjacentReferenceNodeIter->second);
					(leftAdjacentReferenceNodeIter->second)->addOutNode(altNodePtr);
					(rightAdjacentReferenceNodeIter->second)->addInNode(altNodePtr);
				}
			}
		}
	}

	void Graph::addAdditionalNodeEdges()
	{
		Node::SharedPtr nodePtr = this->m_first_node;
		while (nodePtr->getOutNodes().size() > 0)
		{
			for (auto siblingNodePtr : nodePtr->getReferenceOutNode()->getInNodes())
			{
				for (auto rightAdjacentNodePtr : nodePtr->getOutNodes())
				{
					siblingNodePtr->addOutNode(rightAdjacentNodePtr);
					rightAdjacentNodePtr->addInNode(siblingNodePtr);
				}
			}
			nodePtr = nodePtr->getReferenceOutNode();
		}
		nodePtr = this->m_first_node;
		while (nodePtr != nullptr)
		{
			// std::cout << nodePtr->getSequence() << std::endl;
			for (auto allelePtr : nodePtr->getAllelePtrs())
			{
				allelePtr->registerNodePtr(nodePtr);
			}
			nodePtr = nodePtr->getReferenceOutNode();
		}
	}

	Node::SharedPtr Graph::condenseGraph(Node::SharedPtr lastNodePtr)
	{
		Node::SharedPtr nodePtr = lastNodePtr;
		do
		{
			Node::SharedPtr leftAdjacentNodePtr = nodePtr->getReferenceInNode();
			this->m_all_created_nodes.emplace(leftAdjacentNodePtr);
			if (leftAdjacentNodePtr == nullptr)
			{
				std::cout << "Invalid Graph: condenseGraph, position: " << nodePtr->getPosition() << std::endl;
				exit(EXIT_FAILURE);
			}
			if ((leftAdjacentNodePtr->getOutNodes().size() > 1) || (nodePtr->getInNodes().size() > 1))
			{
				nodePtr = leftAdjacentNodePtr;
			}
			else
			{
				nodePtr = Node::mergeNodes(leftAdjacentNodePtr, nodePtr);
			}
			this->m_all_created_nodes.emplace(nodePtr);
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

	void Graph::removeNodePtr(Node* nodePtr)
	{
		std::unordered_set< Node::SharedPtr > allCreatedNodePtrs;
		for (auto iter : m_all_created_nodes)
		{
			if (iter->getID() != nodePtr->getID())
			{
				allCreatedNodePtrs.emplace(nodePtr);
			}
		}
		m_all_created_nodes = allCreatedNodePtrs;
		std::unordered_map< uint32_t, Node::SharedPtr > nodePtrsMap;
		for (auto iter : m_node_ptrs_map)
		{
			if (iter.first != nodePtr->getID())
			{
				nodePtrsMap.emplace(iter.first, iter.second);
			}
		}
		m_node_ptrs_map = nodePtrsMap;
	}

	void getAllPaths(Node::SharedPtr nodePtr, std::vector< Node::SharedPtr > currentPath, int numberOfSibs, std::vector< std::vector< Node::SharedPtr > >& paths)
	{
		currentPath.emplace_back(nodePtr);
		std::unordered_set< Node::SharedPtr > outNodePtrs = nodePtr->getOutNodes();
		if (outNodePtrs.size() == 0)
		{
			paths.emplace_back(currentPath);
		}
		else
		{
			for (auto nextNodePtr : outNodePtrs)
			{
				getAllPaths(nextNodePtr, currentPath, outNodePtrs.size() - 1, paths);
			}
		}
	}

	std::vector< std::vector< Node::SharedPtr > > Graph::generateAllPaths()
	{
		std::vector< std::vector< Node::SharedPtr > > paths;
		Node::SharedPtr firstNode = this->m_first_node;
		std::vector< Node::SharedPtr > path;
		getAllPaths(firstNode, path, 0, paths);
		return paths;
	}

	Region::SharedPtr Graph::getGraphRegion()
	{
		Region::BASED based;
		int startPosition = -1;
		int endPosition = -1;
		std::string referenceID = "";
		for (auto regionPtr : m_graph_regions)
		{
			if (referenceID.size() == 0)
			{
				referenceID = regionPtr->getReferenceID();
				based = regionPtr->getBased();
			}
			if (startPosition < 0 || startPosition > regionPtr->getStartPosition())
			{
				startPosition = regionPtr->getStartPosition();
			}
			if (endPosition < 0 || endPosition > regionPtr->getEndPosition())
			{
				endPosition = regionPtr->getEndPosition();
			}
		}
		return std::make_shared< Region >(referenceID, startPosition, endPosition, based);
	}

	std::string Graph::getReferenceSequence()
	{
		std::string refSequence = "";
		Node::SharedPtr currentRefNode = this->m_first_node;
		while (currentRefNode != nullptr)
		{
			refSequence += currentRefNode->getSequence();
			currentRefNode = currentRefNode->getReferenceOutNode();
		}
		return refSequence;
	}

	void getAllPathsOfLength(Node::SharedPtr nodePtr, std::string currentPath, std::vector< std::string >& paths, uint32_t len)
	{
		if (nodePtr->getSequence().size() < len)
		{
			currentPath += nodePtr->getSequence().substr(len);
			paths.emplace_back(currentPath);
		}
		else
		{
			currentPath += nodePtr->getSequence();
			for (auto nextNodePtr : nodePtr->getOutNodes())
			{
				getAllPathsOfLength(nextNodePtr, currentPath, paths, (nodePtr->getSequence().size() - len));
			}
		}
	}

	std::vector< std::string > Graph::generateAllPathsFromNodesOfLength(Node::SharedPtr nodePtr)
	{
		Node::SharedPtr maxSizeNodePtr = *std::max_element(nodePtr->getOutNodes().begin(), nodePtr->getOutNodes().end(),[](Node::SharedPtr iterNode1,  Node::SharedPtr iterNode2) {
				return iterNode1->getOriginalSequenceSize() < iterNode2->getOriginalSequenceSize();
			});
		std::vector< std::string > paths;
		for (auto outNodePtr : nodePtr->getOutNodes())
		{
			if (maxSizeNodePtr == outNodePtr)
			{
				paths.emplace_back(outNodePtr->getOriginalSequence());
			}
			else
			{
				getAllPathsOfLength(outNodePtr, "", paths, maxSizeNodePtr->getOriginalSequenceSize());
			}
		}
		return paths;
	}

	bool Graph::isNodePrefixAmbiguous(std::string& nodeSequence, Node* comparatorNode, std::unordered_set< Node::SharedPtr >& nodePtrs)
	{
		for (auto nodePtr : nodePtrs)
		{
			auto testNodeSequence = nodePtr->getSequence();
			if (nodePtr.get() == comparatorNode || testNodeSequence.size() < nodeSequence.size())
			{
				continue;
			}
			else if (testNodeSequence.substr(nodeSequence.size()).compare(nodeSequence) == 0)
			{
				return true;
			}
		}
		return false;
	}

	bool Graph::isNodeSuffixAmbiguous(std::string& nodeSequence, Node* comparatorNode, std::unordered_set< Node::SharedPtr >& nodePtrs)
	{
		for (auto nodePtr : nodePtrs)
		{
			auto testNodeSequence = nodePtr->getSequence();
			if (nodePtr.get() == comparatorNode || testNodeSequence.size() < nodeSequence.size())
			{
				continue;
			}
			else if (testNodeSequence.substr(0, nodeSequence.size()).compare(nodeSequence) == 0)
			{
				return true;
			}
		}
		return false;
	}

	void getFullPathFromNode(std::vector< std::string >& paths, std::string currentPath, Node::SharedPtr currentNodePtr)
	{
		currentPath += currentNodePtr->getSequence();
		auto outNodePtrs = currentNodePtr->getOutNodes();
		if (outNodePtrs.size() == 0)
		{
			paths.emplace_back(currentPath);
		}
		else
		{
			for (auto nextNodePtr : outNodePtrs)
			{
				getFullPathFromNode(paths, currentPath, nextNodePtr);
			}
		}
	}

	std::vector< std::string > Graph::getAllPathsAsStrings()
	{
		std::vector< std::string > paths;
		if (this->m_first_node != nullptr)
		{
			getFullPathFromNode(paths, "", this->m_first_node);
		}
		return paths;
	}

	void Graph::printGraphVisOutput()
	{
		std::unordered_set< Node::SharedPtr > printedNodeList;
		std::unordered_set< Node::SharedPtr > printedNodeEdgeList;
		std::cout << "digraph {" << std::endl;
		std::cout << "  rankdir=LR;" << std::endl;
		std::cout << "  splines=ortho;" << std::endl;
		std::vector< Node::SharedPtr > nodePtrs;
		nodePtrs.emplace_back(this->m_first_node);
		char nextNodeID = 'A';
		std::unordered_map< Node::SharedPtr, char > nodePtrToID;
		while (nodePtrs.size() > 0) // print nodes first
		{
			Node::SharedPtr nodePtr = nodePtrs.back();
			nodePtrs.pop_back();
			if (printedNodeList.find(nodePtr) != printedNodeList.end())
			{
				continue;
			}
			printedNodeList.emplace(nodePtr);
			std::string nodeShape = (nodePtr->getAlleleType() == Node::ALLELE_TYPE::REF) ? "box" : "circle";
			auto iter = nodePtrToID.find(nodePtr);
			if (iter == nodePtrToID.end())
			{
				nodePtrToID[nodePtr] = nextNodeID;
				++nextNodeID;
			}
			char nodeID = (char)nodePtrToID.find(nodePtr)->second;
			std::string label = nodePtr->getSequence() + "\n";
			for (auto allelePtr : nodePtr->getAllelePtrs())
			{
				label += std::to_string(reinterpret_cast<std::uintptr_t>(allelePtr.get())) + " [" + std::to_string(allelePtr->getNodePtrs().size()) + "]\n";
			}
			std::cout << "  " << nodeID << " [label=\"" << label << "\" fontname=Arial shape=" << nodeShape << " ]" << std::endl;
			for (auto outNodePtr : nodePtr->getOutNodes())
			{
				nodePtrs.emplace_back(outNodePtr);
			}
		}
		std::cout << std::endl;
		nodePtrs.emplace_back(this->m_first_node);
		while (nodePtrs.size() > 0)
		{
			Node::SharedPtr nodePtr = nodePtrs.back();
			nodePtrs.pop_back();
			if (printedNodeEdgeList.find(nodePtr) != printedNodeEdgeList.end())
			{
				continue;
			}
			printedNodeEdgeList.emplace(nodePtr);
			char nodeID = (char)nodePtrToID.find(nodePtr)->second;
			std::cout << "  " << nodeID << " -> { ";
			int count = 0;
			for (auto outNodePtr : nodePtr->getOutNodes())
			{
				std::string comma = (count++ < nodePtr->getOutNodes().size() - 1) ? "," : "";
				char outNodeID = (char)nodePtrToID.find(outNodePtr)->second;
				std::cout << outNodeID << comma << " ";
				nodePtrs.emplace_back(outNodePtr);
			}
			std::cout << "};" << std::endl;
		}
		std::cout << "}" << std::endl;
	}

	void Graph::clearResources()
	{
		for (auto nodePtr : m_all_created_nodes)
		{
			for (auto allelePtr : nodePtr->getAllelePtrs())
			{
				allelePtr->clearNodePtrs();
			}
			nodePtr->clearAllelePtrs();
			nodePtr->clearInAndOutNodes();
		}
	}

	Graph::SharedPtr Graph::createCopy()
	{
		auto graphPtr = std::make_shared< Graph >();
		graphPtr->m_fasta_reference_ptr = this->m_fasta_reference_ptr;
		graphPtr->m_variant_ptrs = this->m_variant_ptrs;
		graphPtr->m_graph_regions = this->m_graph_regions;
		graphPtr->m_graph_spacing = this->m_graph_spacing;
		graphPtr->m_score_threshold = this->m_score_threshold;
		graphPtr->m_aligned_read_names = this->m_aligned_read_names;
        graphPtr->m_graph_printer_ptr = this->m_graph_printer_ptr;
		graphPtr->m_node_ptrs_map = this->m_node_ptrs_map;
		graphPtr->m_first_node = this->m_first_node;
		return graphPtr;
	}
}
