#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core2/util/Noncopyable.hpp"

#include <memory>

namespace graphite
{
	class Node;
	class Allele : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
		Allele(const std::string& sequence);
		~Allele();

		std::string getSequence() { return this->m_sequence; }
		void registerNodePtr(std::shared_ptr< Node > nodePtr);
		std::shared_ptr< Node > getNodePtr();

	private:
		std::string m_sequence;
		std::shared_ptr< Node > m_node_ptr;

	};
}

#endif// GRAPHITE_ALLELE_H
