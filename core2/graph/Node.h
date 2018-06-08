#ifndef GRAPHITE_NODE_H
#define GRAPHITE_NODE_H

#include "core2/util/Noncopyable.hpp"
#include "core2/util/Types.h"
#include "core2/allele/Allele.h"

#include <string>
#include <memory>
#include <unordered_set>

namespace graphite
{
	class Node : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Node > SharedPtr;
		enum class ALLELE_TYPE { REF = 0, ALT = 1 };
	    Node(const std::string& sequence, position pos, ALLELE_TYPE alleleType);
		Node(char* sequence, uint32_t length, position pos, ALLELE_TYPE alleleType);
		Node();
        ~Node();

		std::vector< Node::SharedPtr > getInNodes();
		std::vector< Node::SharedPtr > getOutNodes();

		ALLELE_TYPE getAlleleType();
		position getPosition();
		std::string getSequence();
		Node::SharedPtr getReferenceInNode();

		void addInNode(Node::SharedPtr node);
		void addOutNode(Node::SharedPtr node);

		void addOverlappingAllelePtr(Allele::SharedPtr allelePtr);
		std::unordered_set< Allele::SharedPtr > getOverlappingAllelePtrs();
		static Node::SharedPtr mergeNodes(Node::SharedPtr firstNodePtr, Node::SharedPtr secondNodePtr);

	private:
		std::unordered_set< Allele::SharedPtr > m_overlapping_allele_ptr_map;
		std::string m_sequence;
		position m_position;
		ALLELE_TYPE m_allele_type;
		Node::SharedPtr m_in_ref_node;
		std::vector< Node::SharedPtr > m_in_nodes;
		std::vector< Node::SharedPtr > m_out_nodes;
	};
}

#endif// GRAPHITE_NODE_H
