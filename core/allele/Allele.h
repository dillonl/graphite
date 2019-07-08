#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core/util/Noncopyable.hpp"
#include "core/sample/Sample.h"
#include "core/util/Types.h"

#include "api/BamAlignment.h"

#include <memory>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

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
		/* void registerNodePtr(std::shared_ptr< Node > nodePtr); */
		/* std::shared_ptr< Node > getNodePtr(); */
		void incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, int score);
        std::unordered_set< std::string > getScoreCountFromAlleleCountType(const std::string& sampleName, AlleleCountType alleleCountType, bool forwardCount);
		void registerNodePtr(std::shared_ptr< Node > nodePtr) { this->m_node_ptrs.emplace(nodePtr); }
		std::unordered_set< std::shared_ptr< Node > > getNodePtrs() { return this->m_node_ptrs; }
		void clearNodePtrs() { this->m_node_ptrs.clear(); }

	private:
		std::string m_sequence;
		/* std::shared_ptr< Node > m_node_ptr; */
        std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_forward_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
        std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_reverse_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
		std::mutex m_counts_lock;
		std::unordered_set< std::shared_ptr< Node > > m_node_ptrs; // you have to make sure to clear this otherwise you will have hanging shared ptrs
	};
}

#endif// GRAPHITE_ALLELE_H
