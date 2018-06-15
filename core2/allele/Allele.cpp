#include "Allele.h"

#include "core2/graph/Node.h"

namespace graphite
{
	Allele::Allele(const std::string& sequence) :
		m_sequence(sequence)
	{
	}

	Allele::~Allele()
	{
	}

	void Allele::registerNodePtr(Node::SharedPtr nodePtr)
	{
		this->m_node_ptr = nodePtr;
	}

	Node::SharedPtr Allele::getNodePtr()
	{
		return this->m_node_ptr;
	}
}
