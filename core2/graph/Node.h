#ifndef GRAPHITE_NODE_H
#define GRAPHITE_NODE_H

#include "core2/util/Noncopyable.hpp"
#include "core2/util/Types.h"
#include "core2/allele/Allele.h"

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
		Node::SharedPtr getReferenceInNode();
		uint32_t getIdenticalPrefixLength();
		uint32_t getIdenticalSuffixLength();

		void addInNode(Node::SharedPtr node);
		void addOutNode(Node::SharedPtr node);
		void setIdenticalPrefixLength(uint32_t prefixLength);
		void setIdenticalSuffixLength(uint32_t suffixLength);

		void addOverlappingAllelePtr(Allele::SharedPtr allelePtr);
		std::unordered_set< Allele::SharedPtr > getOverlappingAllelePtrs();
		void incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, bool isForwardStrand, int score);
		std::unordered_set< std::string > getScoreCountFromAlleleCountType(const std::string& sampleName, AlleleCountType alleleCountType, bool forwardCount);

		static Node::SharedPtr mergeNodes(Node::SharedPtr firstNodePtr, Node::SharedPtr secondNodePtr);

	private:
		static uint32_t s_id; // the static id counter for all nodes
		uint32_t m_id;
		std::unordered_set< Allele::SharedPtr > m_overlapping_allele_ptr_map;
		std::string m_sequence;
		position m_position;
		ALLELE_TYPE m_allele_type;
		Node::SharedPtr m_in_ref_node;
		Node::SharedPtr m_out_ref_node;
		std::unordered_set< Node::SharedPtr > m_in_nodes;
		std::unordered_set< Node::SharedPtr > m_out_nodes;
		uint32_t m_identical_prefix_length; // the length of sequence that is identical to adjacent sequences
		uint32_t m_identical_suffix_length; // the length of sequence that is identical to adjacent sequences
		std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_forward_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
		std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_reverse_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads

	};
}

#endif// GRAPHITE_NODE_H
