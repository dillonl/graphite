#include "Graph.h"

#include "core/util/Types.h"
#include <deque>

namespace graphite
{
	Graph::Graph(FastaReference::SharedPtr fastaReferencePtr, std::vector< Variant::SharedPtr > variantPtrs, uint32_t graphSpacing, bool printGraph) :
		m_fasta_reference_ptr(fastaReferencePtr),
		m_variant_ptrs(variantPtrs),
		m_graph_spacing(graphSpacing),
		m_score_threshold(70),
		m_graph_printer_ptr(nullptr)
	{
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
		getGraphReference(referenceSequence, referenceRegionPtr);
		Node::SharedPtr firstNodePtr;
		Node::SharedPtr lastNodePtr;
		generateReferenceGraphNode(firstNodePtr, lastNodePtr, referenceSequence, referenceRegionPtr);
		addVariantsToGraph(firstNodePtr);
		firstNodePtr = condenseGraph(lastNodePtr);
		setPrefixAndSuffix(firstNodePtr); // calculate prefix and suffix matching sequences
		this->m_first_node = firstNodePtr;
		// compressLargeNodes();
		setRegionPtrs();
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
		Node::SharedPtr prevNode = nullptr;
		static bool firstTime = false;
		for (auto i = 0; i < referenceSequence.size(); ++i)
		{
			Node::SharedPtr node = std::make_shared< Node >(referenceSequence.c_str() + i, 1, nodePosition, Node::ALLELE_TYPE::REF);
			this->m_all_created_nodes.emplace(node);
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
			position variantPosition = variantPtr->getPosition() - 1;
			position inPosition = variantPosition + variantPtr->getReferenceAllelePtr()->getSequence().size() + 1;
			for (uint32_t i = variantPosition + 1; i < inPosition; ++i)
			{
				auto iter = referenceNodePtrPositionMap.find(i);
				if (iter != referenceNodePtrPositionMap.end())
				{
					iter->second->setAllelePtr(variantPtr->getReferenceAllelePtr());
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
				this->m_all_created_nodes.emplace(altNodePtr);
				altNodePtr->setAllelePtr(altAllelePtr);

				position endPosition = altNodePtr->getPosition() + altNodePtr->getSequence().size();
				if (endPositionNodePtrs.find(endPosition) == endPositionNodePtrs.end())
				{
					endPositionNodePtrs[endPosition];
					endPositionNodePtrs[endPosition].emplace_back(altNodePtr);
				}
				else
				{
					endPositionNodePtrs[endPosition].emplace_back(altNodePtr);
				}

				altNodePtr->addInNode(inReferenceIter->second);
				altNodePtr->addOutNode(outReferenceIter->second);
				(inReferenceIter->second)->addOutNode(altNodePtr);
				(outReferenceIter->second)->addInNode(altNodePtr);

				auto iter = endPositionNodePtrs.find(altNodePtr->getPosition());
				if (iter != endPositionNodePtrs.end()) // we subtract one because we want to connect nodes that bump up against us
				{
					for (auto nodePtrs : iter->second)
					{
						altNodePtr->addInNode(nodePtrs);
						nodePtrs->addOutNode(altNodePtr);
					}
				}
			}

		}
	}
	Node::SharedPtr Graph::condenseGraph(Node::SharedPtr lastNodePtr)
	{
		Node::SharedPtr nodePtr = lastNodePtr;
		do
		{
			Node::SharedPtr refNodePtr = nodePtr->getReferenceInNode();
			this->m_all_created_nodes.emplace(refNodePtr);
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

	void Graph::adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent)
	{
		gssw_sse2_disable();
		int8_t* nt_table = gssw_create_nt_table();
		int8_t* mat = gssw_create_score_matrix(matchValue, mismatchValue);
		gssw_graph* graph = gssw_graph_create(this->m_node_ptrs_map.size());

		unordered_map< uint32_t, gssw_node* > gsswNodePtrsMap;
		Node::SharedPtr nodePtr = m_first_node;
		// std::cout << "node from graph: " << nodePtr->getAllelePtr() << std::endl;
		gssw_node* gsswNode = (gssw_node*)gssw_node_create(nodePtr.get(), nodePtr->getID(), nodePtr->getSequence().c_str(), nt_table, mat);
		gsswNodePtrsMap.emplace(nodePtr->getID(), gsswNode);
		gssw_graph_add_node(graph, gsswNode);
		while (nodePtr != nullptr)
		{
			Node::SharedPtr nextRefNodePtr = nullptr;
			for (auto outNodePtr : nodePtr->getOutNodes())
			{
				if (outNodePtr->getAlleleType() == Node::ALLELE_TYPE::REF)
				{
					nextRefNodePtr = outNodePtr;
				}
				gsswNode = (gssw_node*)gssw_node_create(outNodePtr.get(), outNodePtr->getID(), outNodePtr->getSequence().c_str(), nt_table, mat);
				gssw_graph_add_node(graph, gsswNode);
				gsswNodePtrsMap.emplace(outNodePtr->getID(), gsswNode);
			}
			nodePtr = nextRefNodePtr;
		}

		for (auto iter : this->m_node_ptrs_map)
		{
			Node::SharedPtr nodePtr = iter.second;
			gssw_node* gsswNode = gsswNodePtrsMap[nodePtr->getID()];
			for (auto outNodePtr : nodePtr->getOutNodes())
			{
				auto iter = gsswNodePtrsMap.find(outNodePtr->getID());
				if (iter != gsswNodePtrsMap.end())
				{
					gssw_node* gsswOutNode = iter->second;
					gssw_nodes_add_edge(gsswNode, gsswOutNode);
				}
			}
		}

		gssw_graph_fill(graph, bamAlignmentPtr->QueryBases.c_str(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0, 15, 2, true);
		gssw_graph_mapping* gm = gssw_graph_trace_back (graph, bamAlignmentPtr->QueryBases.c_str(), bamAlignmentPtr->QueryBases.size(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0);
		processTraceback(gm, bamAlignmentPtr, samplePtr, !bamAlignmentPtr->IsReverseStrand(), matchValue, mismatchValue, gapOpenValue, gapExtensionValue, referenceTotalScorePercent);
		gssw_graph_mapping_destroy(gm);

		// note that nodes which are referred to in this graph are destroyed as well
		gssw_graph_destroy(graph);

		free(nt_table);
		free(mat);
	}

	void Graph::processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent)
	{
		std::string alignmentName = bamAlignmentPtr->Name + std::to_string(bamAlignmentPtr->IsFirstMate());
		{
			std::lock_guard< std::mutex > l(m_aligned_read_names_mutex);
			if (this->m_aligned_read_names.find(alignmentName) != this->m_aligned_read_names.end())
			{
				return;
			}
			this->m_aligned_read_names.emplace(alignmentName);
		}
		std::vector< std::tuple< Node*, uint32_t > > nodePtrScoreTuples;
		uint32_t totalScore = 0;
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		uint32_t softclipLength = 0;
		std::string fullCigarString = "";
		uint32_t softclipCount = 0;
		bool hasAlternate = false;
		uint32_t firstNodeAlignmentOverlapSize = 0;
		uint32_t lastNodeAlignmentOverlapSize = 0;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			std::string cigarString = "";
			gssw_node* gsswNode = graphMapping->cigar.elements[i].node;
			Node* nodePtr = (Node*)gsswNode->data;
			hasAlternate |= (nodePtr->getAlleleType() == Node::ALLELE_TYPE::ALT);
			int32_t score = 0;
			uint32_t length = 0;
			uint32_t tmpSofclipLength = 0;
			uint32_t nodeLengthAccumulator = 0;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					nodeLengthAccumulator += nc->cigar->elements[j].length;
					score += (matchValue * nc->cigar->elements[j].length);
					break;
				case 'X':
					nodeLengthAccumulator += nc->cigar->elements[j].length;
					score -= (mismatchValue * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					nodeLengthAccumulator += nc->cigar->elements[j].length;
					score -= gapOpenValue;
					score -= (gapExtensionValue * (nc->cigar->elements[j].length -1));
					break;
				case 'S':
					tmpSofclipLength += nc->cigar->elements[j].length;
					++softclipCount;
				default:
					break;
				}
				length += nc->cigar->elements[j].length;
				fullCigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
			}
			if (i == 0)
			{
				firstNodeAlignmentOverlapSize = nodeLengthAccumulator;
			}
			if (i == graphMapping->cigar.length - 1)
			{
				lastNodeAlignmentOverlapSize = nodeLengthAccumulator;
			}
			score = (score < 0) ? 0 : score; // the floor of the mapping score is 0
			softclipLength += tmpSofclipLength;
			float tmpScore = score;
			uint32_t alleleSWScorePercent = (length > 0) ? (tmpScore / (length - tmpSofclipLength)) * 100 : 0;
			std::tuple< Node*, uint32_t > nodePtrScoreTuple(nodePtr, alleleSWScorePercent);
			nodePtrScoreTuples.emplace_back(nodePtrScoreTuple);
			totalScore += score;
		}

		float totalScorePercent = ((float)(totalScore))/((float)(bamAlignmentPtr->Length - softclipLength)) * 100;
		uint32_t count = 0;
		// static std::mutex lo;
		// std::lock_guard< std::mutex > lock(lo);
		for (auto nodePtrScoreTuple : nodePtrScoreTuples)
		{
			bool isFirstNode = count == 0;
			bool isPenultimateNode = count == (nodePtrScoreTuples.size() - 2);
			bool isLastNode = count == (nodePtrScoreTuples.size() - 1);
			Node* nodePtr = std::get< 0 >(nodePtrScoreTuple);
			uint32_t nodeScore = std::get< 1 >(nodePtrScoreTuple);
			bool isAmbiguous = false;
			if (nodePtrScoreTuples.size() > 1)
			{
				if (isFirstNode) // if nodePtr is a candidate for identifying repeats at the beginning of the graph
				{
					std::string firstNodeSequence = nodePtr->getSequence().substr(nodePtr->getSequence().size() - firstNodeAlignmentOverlapSize);
					Node* nextNodePtr = std::get< 0 >(nodePtrScoreTuples[1]);
					std::string concatSequence = nodePtr->getSequence() + nextNodePtr->getSequence();
					for (auto potentialPathNode : nextNodePtr->getInNodes())
					{
						if (nodePtr == potentialPathNode.get()) { continue; }
						std::string potentialConcatSequence = std::string(firstNodeSequence + nextNodePtr->getSequence(), 0, concatSequence.size());
						if (potentialConcatSequence.compare(concatSequence) == 0)
						{
							isAmbiguous = true;
							break;
						}
					}
				}
				else if (isPenultimateNode) // if nodePtr the node before the last node
				{
					Node* lastNodePtr = std::get< 0 >(nodePtrScoreTuples[nodePtrScoreTuples.size() - 1]);
					std::string lastNodeSequence = lastNodePtr->getSequence().substr(0, lastNodeAlignmentOverlapSize);
					if (lastNodePtr->getInNodes().size() > 1 && lastNodePtr->getAlleleType() == Node::ALLELE_TYPE::REF) // if nodePtr is a ref backbone
					{
						std::string concatSequence = nodePtr->getSequence() + lastNodeSequence;
						for (auto potentialPrevPathNode : lastNodePtr->getInNodes())
						{
							if (potentialPrevPathNode.get() == nodePtr) { continue; }
							std::string potentialPathSequence = std::string(potentialPrevPathNode->getSequence() + lastNodePtr->getSequence(), 0, concatSequence.size());
							/*
							static std::mutex lo;
							lo.lock();
							std::cout << "--------------------" << std::endl;
							std::cout << bamAlignmentPtr->Name << std::endl;
							std::cout << concatSequence << std::endl;
							std::cout << potentialPathSequence << std::endl;
							std::cout << "--------------------" << std::endl;
							lo.unlock();
							*/
							if (concatSequence.compare(potentialPathSequence) == 0)
							{
								isAmbiguous = true;
								break;
							}
						}
					}
				}
				else if (isLastNode && nodePtr->getReferenceOutNode() != nullptr && nodePtr->getReferenceOutNode()->getInNodes().size() > 1)
				{
					for (auto potentialNodePtr : nodePtr->getReferenceOutNode()->getInNodes())
					{
						if (nodePtr == potentialNodePtr.get()) { continue; }
						if (nodePtr->getSequence().size() > potentialNodePtr->getSequence().size())
						{
							for (auto nextPotentialNodePtr : potentialNodePtr->getOutNodes())
							{

								std::string potentialSequence = std::string(potentialNodePtr->getSequence() + nextPotentialNodePtr->getSequence(), 0, nodePtr->getSequence().size());
								if (potentialSequence.compare(nodePtr->getSequence()) == 0)
								{
									isAmbiguous = true;
									break;
								}
							}
						}
					}
				}
			}

			bool printed = false;
			// if ((totalScorePercent == referenceTotalScorePercent && hasAlternate) || isAmbiguous)
			if (isAmbiguous)
			{
				nodePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, isForwardStrand, -1);
			}
			else if (totalScorePercent < m_score_threshold || softclipCount > 1)
			{
				nodePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, isForwardStrand, 0);
			}
			else
			{
				nodePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, isForwardStrand, nodeScore);
				if (m_graph_printer_ptr != nullptr)
				{
					m_graph_printer_ptr->registerTraceback(graphMapping, bamAlignmentPtr, totalScorePercent);
					printed = true;
				}
			}
			if (!printed && m_graph_printer_ptr != nullptr)
			{
				m_graph_printer_ptr->registerUnalignedRead(bamAlignmentPtr, fullCigarString, totalScorePercent);
			}
			++count;
		}
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
}
