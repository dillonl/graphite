#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core/util/Noncopyable.hpp"
#include "core/sample/Sample.h"
#include "core/util/Types.h"
#include "core/alignment/Alignment.h"

#include <memory>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

namespace graphite
{
	class Node;
	class SupportingReadInfo : private Noncopyable
	{
	public:
		typedef std::shared_ptr< SupportingReadInfo > SharedPtr;
	    SupportingReadInfo(const std::string& sampleName, const std::string& readName, const std::string& mateReadName, const std::string& cigarString, int nodeScore, int tracebackScore) :
    		m_sample_name(sampleName), m_read_name(readName), m_mate_read_name(mateReadName), m_cigar_string(cigarString), m_node_score(nodeScore), m_traceback_score(tracebackScore)
		{
		}

			static std::string getHeader(const std::string& token) { return "SampleName" + token + "ReadName" + token + "MateReadName" + token + "CigarString" + token + "NodeScore" + token + "TracebackScore"; }
		std::string toString(const std::string& token) { return this->m_sample_name + token + this->m_read_name + token + this->m_mate_read_name + token + this->m_cigar_string + token + std::to_string(this->m_node_score) + token + std::to_string(this->m_traceback_score); }

	private:
		std::string m_sample_name;
		std::string m_read_name;
		std::string m_mate_read_name;
		std::string m_cigar_string;
		int m_node_score;
		int m_traceback_score;
	};

	class Allele : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
		Allele(const std::string& sequence);
		~Allele();

		std::string getSequence() { return this->m_sequence; }
		/* void registerNodePtr(std::shared_ptr< Node > nodePtr); */
		/* std::shared_ptr< Node > getNodePtr(); */
        void incrementScoreCount(Alignment::SharedPtr, Sample::SharedPtr samplePtr, int score);
        std::unordered_set< std::string > getScoreCountFromAlleleCountType(const std::string& sampleName, AlleleCountType alleleCountType, bool forwardCount);
		void registerNodePtr(std::shared_ptr< Node > nodePtr) { this->m_node_ptrs.emplace(nodePtr); }
		std::unordered_set< std::shared_ptr< Node > > getNodePtrs() { return this->m_node_ptrs; }
		void clearNodePtrs() { this->m_node_ptrs.clear(); }
		void pairAllele(Allele::SharedPtr allelePtr);
		void addSemanticLoci(position pos, const std::string& refSequence, const std::string& altSequence);
		std::unordered_map< position, std::unordered_set< std::string > > getSemanticLocations() { return this->m_semantic_locations; }
		void registerSupportingReadInformation(SupportingReadInfo::SharedPtr supportingReadInfo);
		std::vector< SupportingReadInfo::SharedPtr > getSupportingReadInfoPtrs();

	private:
		std::string m_sequence;

		std::unordered_map< position, std::unordered_set< std::string > > m_semantic_locations;
		std::unordered_set< Allele::SharedPtr > m_paired_allele_ptrs;
		/* std::shared_ptr< Node > m_node_ptr; */
        std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_forward_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
        std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > > m_reverse_counts; // map keyed by read SampleName then they are indexed via the AlleleCountType enum value, then add readName to the the unordered set so we are properly counting the reads
		std::mutex m_counts_lock;
		std::unordered_set< std::shared_ptr< Node > > m_node_ptrs; // you have to make sure to clear this otherwise you will have hanging shared ptrs
		std::mutex m_supporting_read_info_mutex;
		std::vector< SupportingReadInfo::SharedPtr > m_supporting_read_info_list;
	};
}

#endif// GRAPHITE_ALLELE_H
