#include "Node.h"

namespace graphite
{
	uint32_t Node::s_id = 0;
	std::mutex Node::s_id_lock;
	Node::Node(const std::string& sequence, position pos, ALLELE_TYPE alleleType) :
		m_sequence(sequence),
		m_position(pos),
		m_allele_type(alleleType),
		m_in_ref_node(nullptr),
		m_original_sequence("")
	{
		s_id_lock.lock();
		m_id = s_id;
		s_id += 1;
		s_id_lock.unlock();
	}

	Node::Node(char* sequence, uint32_t length, position pos, ALLELE_TYPE alleleType) :
		m_sequence(std::string(sequence, length)),
		m_position(pos),
		m_allele_type(alleleType),
		m_in_ref_node(nullptr),
		m_original_sequence("")
	{
		s_id_lock.lock();
		m_id = s_id;
		s_id += 1;
		s_id_lock.unlock();
	}

	Node::Node()
	{
	}

	Node::~Node()
	{
	}

	std::string Node::getSequence()
	{
		return this->m_sequence;
	}

	uint32_t Node::getID()
	{
		return this->m_id;
	}

	Node::SharedPtr Node::getReferenceInNode()
	{
		return this->m_in_ref_node;
	}

	Node::SharedPtr Node::getReferenceOutNode()
	{
		return this->m_out_ref_node;
	}

	std::unordered_set< Node::SharedPtr > Node::getInNodes()
	{
		return this->m_in_nodes;
	}

	std::unordered_set< Node::SharedPtr > Node::getOutNodes()
	{
		return this->m_out_nodes;
	}

	void Node::setIdenticalPrefixLength(uint32_t prefixLength)
	{
		this->m_identical_prefix_length = prefixLength;
	}

	void Node::setIdenticalSuffixLength(uint32_t suffixLength)
	{
		this->m_identical_suffix_length = suffixLength;
	}

	uint32_t Node::getIdenticalPrefixLength()
	{
		return this->m_identical_prefix_length;
	}

	uint32_t Node::getIdenticalSuffixLength()
	{
		return this->m_identical_suffix_length;
	}

	Node::ALLELE_TYPE Node::getAlleleType()
	{
		return this->m_allele_type;
	}

	position Node::getPosition()
	{
		return m_position;
	}

	void Node::addInNode(std::shared_ptr< Node > node)
	{
		if (node->getAlleleType() == ALLELE_TYPE::REF)
		{
			this->m_in_ref_node = node;
		}
		if (this->m_in_nodes.find(node) == this->m_in_nodes.end()) // don't add duplicates
		{
			this->m_in_nodes.emplace(node);
		}
	}

	void Node::addOutNode(std::shared_ptr< Node > node)
	{
		if (node->getAlleleType() == ALLELE_TYPE::REF)
		{
			this->m_out_ref_node = node;
		}
		if (this->m_out_nodes.find(node) == this->m_out_nodes.end()) // don't add duplicates
		{
			this->m_out_nodes.emplace(node);
		}
	}

	void Node::addOverlappingAllelePtr(Allele::SharedPtr allelePtr)
	{
		this->m_overlapping_allele_ptr_map.emplace(allelePtr);
	}

	std::unordered_set< Allele::SharedPtr > Node::getOverlappingAllelePtrs()
	{
		return this->m_overlapping_allele_ptr_map;
	}

	Node::SharedPtr Node::mergeNodes(Node::SharedPtr firstNodePtr, Node::SharedPtr secondNodePtr)
	{
		Node::SharedPtr nodePtr = std::make_shared< Node >(std::string(firstNodePtr->m_sequence + secondNodePtr->m_sequence), firstNodePtr->m_position, firstNodePtr->m_allele_type);
		nodePtr->m_in_nodes = firstNodePtr->m_in_nodes;

		// make sure to set the in and out ref node pointers
		if (nodePtr->m_allele_type == ALLELE_TYPE::REF)
		{
			if (firstNodePtr->m_in_ref_node != nullptr)
			{
				firstNodePtr->m_in_ref_node->m_out_ref_node = nodePtr;
			}
			if (secondNodePtr->m_out_ref_node != nullptr)
			{
				firstNodePtr->m_out_ref_node->m_in_ref_node = nodePtr;
			}
		}

		// redirect the in_nodes of firstNodePtr to point to nodePtr
		for (auto inNodePtr : nodePtr->m_in_nodes)
		{
			inNodePtr->m_out_nodes.erase(firstNodePtr);
			inNodePtr->m_out_nodes.emplace(nodePtr);
		}

		nodePtr->m_out_nodes = secondNodePtr->m_out_nodes;
		for (auto outNodePtr : nodePtr->m_out_nodes)
		{
			outNodePtr->m_in_nodes.erase(secondNodePtr);
			outNodePtr->m_in_nodes.emplace(nodePtr);
		}
		nodePtr->m_overlapping_allele_ptr_map = firstNodePtr->m_overlapping_allele_ptr_map;
		nodePtr->m_in_ref_node = firstNodePtr->m_in_ref_node;
		nodePtr->m_out_ref_node = secondNodePtr->m_out_ref_node;
		for (auto allelePtr : secondNodePtr->m_overlapping_allele_ptr_map)
		{
			nodePtr->m_overlapping_allele_ptr_map.emplace(allelePtr);
		}
		if (firstNodePtr->getAllelePtrs().size() > 0)
		{
			for (auto allelePtr : firstNodePtr->getAllelePtrs())
			{
				nodePtr->registerAllelePtr(allelePtr);
			}
		}
		else if (secondNodePtr->getAllelePtrs().size() > 0)
		{
			for (auto allelePtr : secondNodePtr->getAllelePtrs())
			{
				nodePtr->registerAllelePtr(allelePtr);
			}
		}
		return nodePtr;
	}

	/*
	void Node::incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, int score)
	{
		if (this->m_allele_ptr != nullptr)
		{
			this->m_allele_ptr->incrementScoreCount(bamAlignmentPtr, samplePtr, isForwardStrand, score);
		}
	}
	*/

	void Node::clearInAndOutNodes()
	{
		this->m_in_ref_node = nullptr;
		this->m_out_ref_node = nullptr;
		this->m_out_nodes.clear();
		this->m_in_nodes.clear();
	}

	void Node::setCompressedSequence(const std::string& sequence)
	{
		this->m_original_sequence = this->m_sequence;
		this->m_sequence = sequence;
	}

	std::string Node::getOriginalSequence()
	{
		if (this->m_original_sequence.size() > 0)
		{
			return this->m_original_sequence;
		}
		else
		{
			return this->m_sequence;
		}
	}

	uint32_t Node::getOriginalSequenceSize()
	{
		if (this->m_original_sequence.size() > 0)
		{
			return this->m_original_sequence.size();
		}
		else
		{
			return this->m_sequence.size();
		}
	}

	void Node::registerAllelePtr(Allele::SharedPtr allelePtr)
	{
		this->m_allele_ptrs.emplace(allelePtr);
	}


}
