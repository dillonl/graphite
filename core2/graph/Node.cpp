#include "Node.h"

namespace graphite
{
	Node::Node(const std::string& sequence, position pos, ALLELE_TYPE alleleType) :
		m_sequence(sequence),
		m_position(pos),
		m_allele_type(alleleType),
		m_in_ref_node(nullptr)
	{
	}

	Node::Node(char* sequence, uint32_t length, position pos, ALLELE_TYPE alleleType) :
		m_sequence(std::string(sequence, length)),
		m_position(pos),
		m_allele_type(alleleType),
		m_in_ref_node(nullptr)
	{
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

	Node::SharedPtr Node::getReferenceInNode()
	{
		return this->m_in_ref_node;
	}

	std::vector< Node::SharedPtr > Node::getInNodes()
	{
		return this->m_in_nodes;
	}

	std::vector< Node::SharedPtr > Node::getOutNodes()
	{
		return this->m_out_nodes;
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
		this->m_in_nodes.emplace_back(node);
	}

	void Node::addOutNode(std::shared_ptr< Node > node)
	{
		this->m_out_nodes.emplace_back(node);
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
		nodePtr->m_in_nodes = secondNodePtr->m_out_nodes;
		nodePtr->m_overlapping_allele_ptr_map = firstNodePtr->m_overlapping_allele_ptr_map;
		for (auto allelePtr : secondNodePtr->m_overlapping_allele_ptr_map)
		{
			nodePtr->m_overlapping_allele_ptr_map.emplace(allelePtr);
		}
		return nodePtr;
	}
}
