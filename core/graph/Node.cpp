#include "Node.h"

namespace graphite
{
	uint32_t Node::s_id = 0;
	Node::Node(const std::string& sequence, position pos, ALLELE_TYPE alleleType) :
		m_sequence(sequence),
		m_position(pos),
		m_allele_type(alleleType),
		m_in_ref_node(nullptr),
		m_identical_prefix_length(0),
		m_identical_suffix_length(0),
		m_id(s_id),
		m_original_sequence("")
	{
		s_id += 1;
	}

	Node::Node(char* sequence, uint32_t length, position pos, ALLELE_TYPE alleleType) :
		m_sequence(std::string(sequence, length)),
		m_position(pos),
		m_allele_type(alleleType),
		m_in_ref_node(nullptr),
		m_identical_prefix_length(0),
		m_identical_suffix_length(0),
		m_id(s_id),
		m_original_sequence("")
	{
		s_id += 1;
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
		this->m_in_nodes.emplace(node);
	}

	void Node::addOutNode(std::shared_ptr< Node > node)
	{
		if (node->getAlleleType() == ALLELE_TYPE::REF)
		{
			this->m_out_ref_node = node;
		}
		this->m_out_nodes.emplace(node);
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
		return nodePtr;
	}

	void Node::incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, int score)
	{
		bool nodeSkipped = true;
		for (auto allelePtr : this->m_overlapping_allele_ptr_map)
		{
			nodeSkipped = false;
			allelePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, isForwardStrand, score);
		}
		/*
		std::string nodeType = (this->getAlleleType() == Node::ALLELE_TYPE::REF) ? "REF" : "ALT";
		if (nodeSkipped)
		{
			// std::cout << "node skipped: " << this->m_id << std::endl;
			std::cout << "node skipped: " << nodeType << std::endl;
		}
		else
		{
			// std::cout << nodeType << " " << score << " " << AlleleCountTypeToString(scoreToAlleleCountType(score)) << std::endl;
			std::cout << "node not skipped: " << nodeType << std::endl;
		}
		*/
	}

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


}
