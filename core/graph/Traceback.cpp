#include "Traceback.h"

namespace graphite
{
	Traceback::Traceback()
	{
	}

	Traceback::~Traceback()
	{
	}

	void Traceback::processTraceback(gssw_graph_mapping* graphMapping, Alignment::SharedPtr alignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue)
	{
		this->m_total_score = 0;
		this->m_traceback_nodes.clear();
		this->m_number_of_softclips = 0;
		uint32_t totalSoftclipLength = 0;
		int32_t totalScore = 0;
		std::string graphCigarString = "";
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
				graphCigarString += nc->cigar->elements[j].type + std::to_string(nc->cigar->elements[j].length);
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
		if (alignmentPtr->getLength() > 0)
		{
			this->m_total_score = ((float)totalScore / (float)((alignmentPtr->getLength() - totalSoftclipLength) * matchValue)) * 100;
		}
		if (this->m_total_score >= 80 && this->m_number_of_softclips <= 1 && totalSoftclipLength < (alignmentPtr->getLength() * 0.3))
		{
			this->incrementAlleleCounts(alignmentPtr, samplePtr, graphCigarString);
		}
	}

	void Traceback::incrementAlleleCounts(Alignment::SharedPtr alignmentPtr, Sample::SharedPtr samplePtr, const std::string& graphCigarString)
	{
		for (auto tracebackNodePtr : this->m_traceback_nodes)
		{
			if (tracebackNodePtr->getScore() >= 70 && !nodeHasFlankingMismatches(tracebackNodePtr, alignmentPtr))
			{
				auto nodeScore = (isNodeSequenceAmbiguous(tracebackNodePtr)) ? -1 : tracebackNodePtr->getScore();
				// auto nodeScore = tracebackNodePtr->getScore();
				for (auto nodeAllelePtr : tracebackNodePtr->getNodePtr()->getAllelePtrs())
				{
					if (fullAlleleInTraceback(nodeAllelePtr, tracebackNodePtr->getNodePtr()))
					{
						nodeAllelePtr->incrementScoreCount(alignmentPtr, nodeScore);
						if (nodeScore != 0)
						{
							std::string alignmentName = alignmentPtr->getReadName() +	std::to_string(!alignmentPtr->getIsFirstMate() + 1);
							std::string mateAlignmentName = alignmentPtr->getReadName() +	std::to_string((alignmentPtr->getIsFirstMate()) + 1);
							auto supportingReadInfo = std::make_shared< SupportingReadInfo >(samplePtr->getName(), alignmentName, mateAlignmentName, graphCigarString, nodeScore, this->m_total_score);
						    nodeAllelePtr->registerSupportingReadInformation(supportingReadInfo);
						}
					}
				}
			}
		}
	}

	bool Traceback::nodeHasFlankingMismatches(TracebackNode::SharedPtr tracebackNodePtr, Alignment::SharedPtr alignmentPtr)
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
		bool recordedRef = false;
		auto prevTracebackNodePtr = tracebackNodePtr->getPrevTracebackNodePtr();
		auto nextTracebackNodePtr = tracebackNodePtr->getNextTracebackNodePtr();

		if (prevTracebackNodePtr == nullptr && nextTracebackNodePtr != nullptr) // I am the first node
		{
			uint32_t matchCount = tracebackNodePtr->getMatchCountFromEnd();
			std::string matchSequence = tracebackNodePtr->getNodePtr()->getSequence().substr(tracebackNodePtr->getNodePtr()->getSequence().size() - matchCount);
			for (auto siblingNodePtr : nextTracebackNodePtr->getNodePtr()->getInNodes())
			{
				if (tracebackNodePtr->getNodePtr() == siblingNodePtr.get() || siblingNodePtr->getSequence().size() < matchCount) { continue; } // we don't want to compare ourselves!
				std::string siblingMatchSequence = siblingNodePtr->getSequence().substr(tracebackNodePtr->getNodePtr()->getSequence().size() - matchCount);
				if (siblingMatchSequence.compare(matchSequence) == 0)
				{
					return true;
				}
			}
		}
		if (nextTracebackNodePtr == nullptr && prevTracebackNodePtr != nullptr) // I am the last node
		{
			uint32_t matchCount = tracebackNodePtr->getMatchCountFromStart();
			std::string matchSequence = tracebackNodePtr->getNodePtr()->getSequence().substr(0, matchCount);
			for (auto siblingNodePtr : prevTracebackNodePtr->getNodePtr()->getOutNodes())
			{
				if (tracebackNodePtr->getNodePtr() == siblingNodePtr.get() || siblingNodePtr->getSequence().size() < matchCount) { continue; } // we don't want to compare ourselves!
				std::string siblingMatchSequence = siblingNodePtr->getSequence().substr(0, matchCount);
				if (siblingMatchSequence.compare(matchSequence) == 0)
				{
					return true;
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
