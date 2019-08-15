#include "Traceback.h"

namespace graphite
{
	Traceback::Traceback(gssw_graph* graph, gssw_graph_mapping* gm, int8_t* nt_table, int8_t* mat, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent) :
		m_graph(graph),
		m_gm(gm),
		m_nt_table(nt_table),
		m_mat(mat),
		m_bam_alignment_ptr(bamAlignmentPtr),
		m_sample_ptr(samplePtr),
		m_match_value(matchValue),
		m_mismatch_value(mismatchValue),
		m_gap_open_value(gapOpenValue),
		m_gap_extension_value(gapExtensionValue),
		m_reference_total_score_percent(referenceTotalScorePercent)
	{
	}

	Traceback::~Traceback()
	{
		gssw_graph_mapping_destroy(m_gm);

		// note that nodes which are referred to in this graph are destroyed as well
		gssw_graph_destroy(m_graph);

		free(m_nt_table);
		free(m_mat);
	}

	void Traceback::processTraceback()
	{
		gssw_graph_fill(m_graph, m_bam_alignment_ptr->QueryBases.c_str(), m_nt_table, m_mat, m_gap_open_value, m_gap_extension_value, 0, 0, 15, 2, true);
		this->m_gm = gssw_graph_trace_back(m_graph, m_bam_alignment_ptr->QueryBases.c_str(), m_bam_alignment_ptr->QueryBases.size(), m_nt_table, m_mat, m_gap_open_value, m_gap_extension_value, 0, 0);

		this->m_total_score = 0;
		this->m_traceback_nodes.clear();
		this->m_number_of_softclips = 0;
		uint32_t totalSoftclipLength = 0;
		int32_t totalScore = 0;
		gssw_node_cigar* nc = this->m_gm->cigar.elements;
		TracebackNode::SharedPtr prevTracebackNodePtr = nullptr;
		for (int i = 0; i < this->m_gm->cigar.length; ++i, ++nc)
		{
			auto tracebackNodePtr = std::make_shared< TracebackNode >();
			gssw_node* gsswNode = this->m_gm->cigar.elements[i].node;
			Node* nodePtr = (Node*)gsswNode->data;
			int32_t nodeScore = 0;
			uint32_t nodeLength = 0;
			uint32_t nodeSoftclipLength = 0;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				std::tuple< uint32_t, char > cigarComponent = std::make_tuple(nc->cigar->elements[j].length, nc->cigar->elements[j].type);
				tracebackNodePtr->addCigarComponent(cigarComponent);
				nodeLength += nc->cigar->elements[j].length;
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					nodeScore += (this->m_match_value * nc->cigar->elements[j].length);
					break;
				case 'X':
					nodeScore -= (m_mismatch_value * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					nodeScore -= m_gap_open_value;
					nodeScore -= (m_gap_extension_value * (nc->cigar->elements[j].length -1));
					break;
				case 'S':
					nodeSoftclipLength += nc->cigar->elements[j].length;
					++this->m_number_of_softclips;
				default:
					break;
				}
			}
			nodeScore = (nodeScore < 0) ? 0 : nodeScore; // the floor of the mapping score is 0
			totalSoftclipLength += nodeSoftclipLength;
			int32_t nodeScorePercent = (nodeLength > 0) ? ((float)nodeScore / ((float)(nodeLength - nodeSoftclipLength) * m_match_value)) * 100 : 0;
			totalScore += nodeScore;
			tracebackNodePtr->setNodePtr(nodePtr);
			tracebackNodePtr->setPrevTracebackNodePtr(prevTracebackNodePtr);
			if (prevTracebackNodePtr != nullptr)
			{
				prevTracebackNodePtr->setNextTracebackNodePtr(tracebackNodePtr);
			}
			tracebackNodePtr->setNodeScore(nodeScorePercent);
			prevTracebackNodePtr = tracebackNodePtr;
			tracebackNodePtr->setNextTracebackNodePtr(nullptr); // this will get set on the next time around unless it's the last node, then we want it nullptr
			this->m_traceback_nodes.emplace_back(tracebackNodePtr);
		}
		if (m_bam_alignment_ptr->QueryBases.size() > 0)
		{
			this->m_total_score = ((float)totalScore / (float)((this->m_bam_alignment_ptr->QueryBases.size() - totalSoftclipLength) * m_match_value)) * 100;
		}
		if (this->m_total_score >= 90 && this->m_number_of_softclips <= 1 && totalSoftclipLength < (this->m_bam_alignment_ptr->QueryBases.size() * 0.3))
		{
			this->incrementAlleleCounts(this->m_bam_alignment_ptr, this->m_sample_ptr);
		}
	}

	void Traceback::incrementAlleleCounts(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr)
	{
		for (auto tracebackNodePtr : this->m_traceback_nodes)
		{
			if (tracebackNodePtr->getScore() >= 70)// && !nodeHasFlankingMismatches(tracebackNodePtr, bamAlignmentPtr))
			{
				auto nodeScore = (isNodeSequenceAmbiguous(tracebackNodePtr)) ? -1 : tracebackNodePtr->getScore();
				for (auto nodeAllelePtr : tracebackNodePtr->getNodePtr()->getAllelePtrs())
				{
					if (fullAlleleInTraceback(nodeAllelePtr, tracebackNodePtr->getNodePtr()))
					{
						nodeAllelePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, nodeScore);
					}
				}
			}
		}
	}

	bool Traceback::nodeHasFlankingMismatches(TracebackNode::SharedPtr tracebackNodePtr, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr)
	{
		auto prevTracebackNodePtr = tracebackNodePtr->getPrevTracebackNodePtr();
		auto nextTracebackNodePtr = tracebackNodePtr->getNextTracebackNodePtr();
		if (prevTracebackNodePtr != nullptr)
		{
			auto cigar = prevTracebackNodePtr->getCigar();
			if (cigar.size() > 0 && std::get< 1 >(cigar[cigar.size() - 1]) != 'M')
			{
				return true;
			}
		}
		if (nextTracebackNodePtr != nullptr)
		{
			auto cigar = nextTracebackNodePtr->getCigar();
			if (cigar.size() > 0 && std::get< 1 >(cigar[0]) != 'M')
			{
				return true;
			}
		}
		return false;
	}

	bool Traceback::isNodeSequenceAmbiguous(TracebackNode::SharedPtr tracebackNodePtr)
	{
		auto prevTracebackNodePtr = tracebackNodePtr->getPrevTracebackNodePtr();
		auto nextTracebackNodePtr = tracebackNodePtr->getNextTracebackNodePtr();
		auto matchCount = tracebackNodePtr->getMatchCount();
		if (prevTracebackNodePtr == nullptr && nextTracebackNodePtr != nullptr)
		{
			std::string nodeSequence = tracebackNodePtr->getNodePtr()->getSequence().substr(tracebackNodePtr->getNodePtr()->getSequence().size() - matchCount);
			for (auto siblingNodePtr : nextTracebackNodePtr->getNodePtr()->getInNodes())
			{
				if (tracebackNodePtr->getNodePtr() == siblingNodePtr.get()) { continue; }
				if (siblingNodePtr->getSequence().size() >= matchCount)
				{
					std::string siblingNodeSequence = siblingNodePtr->getSequence().substr(siblingNodePtr->getSequence().size() - matchCount);
					if (nodeSequence.compare(siblingNodeSequence) == 0)
					{
						return true;
					}
				}
			}
		}
	    else if (nextTracebackNodePtr == nullptr && prevTracebackNodePtr != nullptr)
		{
			std::string nodeSequence = tracebackNodePtr->getNodePtr()->getSequence().substr(0, matchCount);
			for (auto siblingNodePtr : prevTracebackNodePtr->getNodePtr()->getOutNodes())
			{
				if (tracebackNodePtr->getNodePtr() == siblingNodePtr.get()) { continue; }
				if (siblingNodePtr->getSequence().size() >= matchCount)
				{
					std::string siblingNodeSequence = siblingNodePtr->getSequence().substr(0, siblingNodePtr->getSequence().size());
					if (nodeSequence.compare(siblingNodeSequence) == 0)
					{
						return true;
					}
				}
			}
		}
		return false;
	}
	bool Traceback::fullAlleleInTraceback(Allele::SharedPtr allelePtr, Node* nodePtr)
	{
		std::string nodeSequence = nodePtr->getSequence();
		if (nodePtr->getAlleleType() == Node::ALLELE_TYPE::REF // alt nodes do not split up alleles so we don't have to check
			&& allelePtr->getSequence().compare(nodeSequence) != 0) // if the allele and node have different sequences then we need to check
		{
			std::string piecemealAlleleSequence = "";
			for (auto tracebackNodePtr : this->m_traceback_nodes)
			{
				auto tracebackRawNodePtr = tracebackNodePtr->getNodePtr();
				if (tracebackRawNodePtr->hasAllelePtr(allelePtr))
				{
					piecemealAlleleSequence += tracebackRawNodePtr->getSequence();
				}
				else if (piecemealAlleleSequence.size() > 0) // if we have started building the piecemealAlleleSequence and stopped matching then we are done building the sequence and we can break out, this is just a heuristic
				{
					break;
				}
			}
			bool isFullAllelePresent = (nodeSequence.compare(allelePtr->getSequence()) == 0);
			return isFullAllelePresent; // after we have constructed the whole sequence then they should be equal if the traceback goes through the entire allele
		}
		return true;
	}

	void Traceback::printTraceback()
	{
		for (auto tracebackNodePtr : m_traceback_nodes)
		{
			std::cout << tracebackNodePtr->getNodePtr()->getSequence() << "|";
		}
		std::cout << std::endl;
	}
}
