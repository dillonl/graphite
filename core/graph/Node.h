#ifndef GRAPHITE_NODE_H
#define GRAPHITE_NODE_H

#include "core/util/Noncopyable.hpp"
#include "core/util/Types.h"
#include "core/allele/Allele.h"
#include "core/sample/Sample.h"

#include "api/BamAlignment.h"

#include <string>
#include <memory>
#include <unordered_set>
#include <unordered_map>

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

		std::unordered_set< Node::SharedPtr > getInNodes();
		std::unordered_set< Node::SharedPtr > getOutNodes();
		uint32_t getID();

		ALLELE_TYPE getAlleleType();
		position getPosition();
		std::string getSequence();
		void setCompressedSequence(const std::string& sequence);
		Node::SharedPtr getReferenceInNode();
		Node::SharedPtr getReferenceOutNode();
		uint32_t getIdenticalPrefixLength();
		uint32_t getIdenticalSuffixLength();

		void addInNode(Node::SharedPtr node);
		void addOutNode(Node::SharedPtr node);
		void setIdenticalPrefixLength(uint32_t prefixLength);
		void setIdenticalSuffixLength(uint32_t suffixLength);

		void addOverlappingAllelePtr(Allele::SharedPtr allelePtr);
		std::unordered_set< Allele::SharedPtr > getOverlappingAllelePtrs();
		void incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, int score);
		std::string getOriginalSequence();
		uint32_t getOriginalSequenceSize();

		static Node::SharedPtr mergeNodes(Node::SharedPtr firstNodePtr, Node::SharedPtr secondNodePtr);

		void clearInAndOutNodes();

	private:
		static uint32_t s_id; // the static id counter for all nodes
		uint32_t m_id;
		std::unordered_set< Allele::SharedPtr > m_overlapping_allele_ptr_map;
		std::string m_sequence;
		std::string m_original_sequence;
		position m_position;
		ALLELE_TYPE m_allele_type;
		Node::SharedPtr m_in_ref_node;
		Node::SharedPtr m_out_ref_node;
		std::unordered_set< Node::SharedPtr > m_in_nodes;
		std::unordered_set< Node::SharedPtr > m_out_nodes;
		uint32_t m_identical_prefix_length; // the length of sequence that is identical to adjacent sequences
		uint32_t m_identical_suffix_length; // the length of sequence that is identical to adjacent sequences

	};
}

#endif// GRAPHITE_NODE_H
