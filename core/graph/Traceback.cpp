#include "Traceback.h"

namespace graphite
{
	Traceback::Traceback()
	{
	}

	Traceback::~Traceback()
	{
	}

	void Traceback::processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue, float referenceTotalScorePercent)
	{
		this->m_total_score = 0;
		this->m_traceback_nodes.clear();
		this->m_number_of_softclips = 0;
		uint32_t totalSoftclipLength = 0;
		int32_t totalScore = 0;
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		TracebackNode::SharedPtr prevTracebackNodePtr = nullptr;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			auto tracebackNodePtr = std::make_shared< TracebackNode >();
			gssw_node* gsswNode = graphMapping->cigar.elements[i].node;
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
					nodeScore += (matchValue * nc->cigar->elements[j].length);
					break;
				case 'X':
					nodeScore -= (mismatchValue * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					nodeScore -= gapOpenValue;
					nodeScore -= (gapExtensionValue * (nc->cigar->elements[j].length -1));
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
			// int32_t nodeScorePercent = ((float)nodeScore / (float)((nodeLength - nodeSoftclipLength) * matchValue)) * 100;
			int32_t nodeScorePercent = (nodeLength > 0) ? ((float)nodeScore / ((float)(nodeLength - nodeSoftclipLength) * matchValue)) * 100 : 0;
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
		if (bamAlignmentPtr->QueryBases.size() > 0)
		{
			this->m_total_score = ((float)totalScore / (float)((bamAlignmentPtr->QueryBases.size() - totalSoftclipLength) * matchValue)) * 100;
		}
		if (this->m_total_score >= 90 && this->m_number_of_softclips <= 1 && totalSoftclipLength < (bamAlignmentPtr->QueryBases.size() * 0.1) && this->m_traceback_nodes.size() > 1)
		{
			this->incrementAlleleCounts(bamAlignmentPtr, samplePtr);
		}
		else
		{
			/*
			std::lock_guard< std::mutex > l(s_traceback_mutex);
			std::cout << "----------" << std::endl;
			std::cout << bamAlignmentPtr->QueryBases << std::endl;
			std::cout << "failed score and/or softclips: score: " << this->m_total_score << " sc: " << this->m_number_of_softclips << " scl: " << totalSoftclipLength << std::endl;
			std::cout << "----------" << std::endl;
			*/
		}
	}

	void Traceback::incrementAlleleCounts(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr)
	{
		for (auto tracebackNodePtr : this->m_traceback_nodes)
		{
			if (tracebackNodePtr->getScore() >= 90 &&
				!nodeHasFlankingMismatches(tracebackNodePtr, bamAlignmentPtr) &&
				!isNodeSequenceAmbiguous(tracebackNodePtr))
			{
				for (auto nodeAllelePtr : tracebackNodePtr->getNodePtr()->getAllelePtrs())
				{
					if (fullAlleleInTraceback(nodeAllelePtr, tracebackNodePtr->getNodePtr()))
					{
						nodeAllelePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, tracebackNodePtr->getScore());
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
				/*
				std::lock_guard< std::mutex > l(s_traceback_mutex);
				std::cout << "----------" << std::endl;
				std::cout << bamAlignmentPtr->QueryBases << std::endl;
				std::cout << "has left flanking mismatches" << std::endl;
				std::cout << prevTracebackNodePtr->printCigarString() << std::endl;
				std::cout << "----------" << std::endl;
				*/
				return true;
			}
		}
		if (nextTracebackNodePtr != nullptr)
		{
			auto cigar = nextTracebackNodePtr->getCigar();
			if (cigar.size() > 0 && std::get< 1 >(cigar[0]) != 'M')
			{
				/*
				std::lock_guard< std::mutex > l(s_traceback_mutex);
				std::cout << "----------" << std::endl;
				std::cout << bamAlignmentPtr->QueryBases << std::endl;
				std::cout << "has right flanking mismatches" << std::endl;
				std::cout << nextTracebackNodePtr->printCigarString() << std::endl;
				std::cout << "----------" << std::endl;
				*/
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
						/*
						std::lock_guard< std::mutex > l(s_traceback_mutex);
						std::cout << "----------" << std::endl;
						std::cout << "node sequence: " << nodeSequence << std::endl;
						std::cout << "sib node seq : " << siblingNodePtr->getSequence() << std::endl;
						std::cout << "beginning node sequence is ambiguous" << std::endl;
						std::cout << "----------" << std::endl;
						*/
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
						/*
						std::lock_guard< std::mutex > l(s_traceback_mutex);
						std::cout << "----------" << std::endl;
						std::cout << "ending node sequence is ambiguous" << std::endl;
						std::cout << "----------" << std::endl;
						*/
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
			if (!isFullAllelePresent)
			{
				/*
				std::lock_guard< std::mutex > l(s_traceback_mutex);
				std::cout << "----------" << std::endl;
				std::cout << "full allele not present" << std::endl;
				std::cout << "----------" << std::endl;
				*/
			}
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
