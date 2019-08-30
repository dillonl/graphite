#pragma once

#include "core/util/Noncopyable.hpp"
#include "Node.h"

#include "api/BamAlignment.h"
#include "gssw.h"

#include <memory>
#include <unordered_map>

/*
  What are all the things we want to check:
  1. total score is high enough
  2. no more than 1 softclip region

  3. individual node score is high enough
  4. make sure traceback travels through all nodes in the allele being considered (if it's at the beginning or end of traceback maybe mark as ambiguous)
  5. make sure flanking nodes don't have mismatches at the edges
  6. if first or last node check for ambiguousness
*/

namespace graphite
{
	class TracebackNode : private Noncopyable
	{
	public:
		typedef std::shared_ptr< TracebackNode > SharedPtr;
	    TracebackNode() {}
		~TracebackNode() {}

		void setPrevTracebackNodePtr(TracebackNode::SharedPtr tracebackNodePtr) { this->m_prev_tracebacknode_ptr = tracebackNodePtr; }
		void setNextTracebackNodePtr(TracebackNode::SharedPtr tracebackNodePtr) { this->m_next_tracebacknode_ptr = tracebackNodePtr; }
		void setNodePtr(Node* nodePtr) { this->m_node_ptr = nodePtr; }
		void setNodeScore(int32_t score) {this->m_score = score;}
		void addCigarComponent(std::tuple< uint32_t, char > cigarComponent) { this->m_cigar.emplace_back(cigarComponent); }
		uint32_t printCigarString()
		{
			for (auto iter : this->m_cigar)
			{
				std::cout << std::to_string(std::get< 0 >(iter)) << std::get< 1 >(iter);
			}
			std::cout << std::endl;
		}

		uint32_t getMatchCountFromStart()
		{
			uint32_t count = 0;
			for (auto iter : this->m_cigar)
			{
				char cigarElem = std::get< 1 >(iter);
				if (cigarElem != 'M') // once we stop matching then break
				{
					break;
				}
				count += std::get< 0 >(iter);
			}
			return count;
		}

		uint32_t getMatchCountFromEnd()
		{
			uint32_t count = 0;
			for (auto iter = this->m_cigar.rbegin(); iter != this->m_cigar.rend(); ++iter)
			{
				char cigarElem = std::get< 1 >(*iter);
				if (cigarElem != 'M') // once we stop matching then break
				{
					break;
				}
				count += std::get< 0 >(*iter);
			}
			return count;
		}

		Node* getNodePtr() { return this->m_node_ptr; }
		TracebackNode::SharedPtr getPrevTracebackNodePtr() { return this->m_prev_tracebacknode_ptr; }
		TracebackNode::SharedPtr getNextTracebackNodePtr() { return this->m_next_tracebacknode_ptr; }
		int32_t getScore() { return this->m_score; }
		std::vector< std::tuple< uint32_t, char > > getCigar() { return this->m_cigar; }

	private:
		Node* m_node_ptr;
		TracebackNode::SharedPtr m_prev_tracebacknode_ptr;
		TracebackNode::SharedPtr m_next_tracebacknode_ptr;
		int32_t m_score;
		std::vector< std::tuple< uint32_t, char > > m_cigar;
	};

    class Traceback : private Noncopyable
	{
	public:
        typedef std::shared_ptr< Traceback > SharedPtr;
		Traceback();
		~Traceback();

		void processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent);
		void printTraceback();

	private:
		void incrementAlleleCounts(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr);
		bool nodeHasFlankingMismatches(TracebackNode::SharedPtr tracebackNodePtr, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr);
		bool isNodeSequenceAmbiguous(TracebackNode::SharedPtr tracebackNodePtr);
		bool fullAlleleInTraceback(Allele::SharedPtr allelePtr, Node* nodePtr);
		std::shared_ptr< BamTools::BamAlignment > m_bam_alignment_ptr;
		Sample::SharedPtr m_sample_ptr;
		int32_t m_total_score;
		uint32_t m_number_of_softclips; // this is not a count of the bases but rather the number of times a sc occurres in the traceback
		std::vector< TracebackNode::SharedPtr > m_traceback_nodes;
	};
}
