#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core2/util/Noncopyable.hpp"
#include "core2/sample/Sample.h"
#include "core2/util/Types.h"

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
		void incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, int score);
        std::unordered_set< std::string > getScoreCountFromAlleleCountType(const std::string& sampleName, AlleleCountType alleleCountType, bool forwardCount);

	private:
		std::string m_sequence;
		/* std::shared_ptr< Node > m_node_ptr; */
        std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_forward_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
        std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_reverse_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
		std::mutex m_counts_lock;
	};
}

#endif// GRAPHITE_ALLELE_H
