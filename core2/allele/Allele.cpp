#include "Allele.h"

#include "core2/graph/Node.h"

namespace graphite
{
	Allele::Allele(const std::string& sequence) :
		m_sequence(sequence)
	{
	}

	Allele::~Allele()
	{
	}

	/*
	void Allele::registerNodePtr(Node::SharedPtr nodePtr)
	{
		this->m_node_ptr = nodePtr;
	}
	*/

	/*
	Node::SharedPtr Allele::getNodePtr()
	{
		return this->m_node_ptr;
	}
	*/

	void Allele::incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, int score)
	{
		size_t alleleCountType = (size_t)scoreToAlleleCountType(score);
		if (score < 0)
		{
			alleleCountType = (size_t)AlleleCountType::Ambiguous;
		}
		std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > >* counts = (isForwardStrand) ? &this->m_forward_counts : &this->m_reverse_counts;

		// std::string nodeType = (this->getAlleleType() == Node::ALLELE_TYPE::REF) ? "REF" : "ALT";
		// std::cout << nodeType << " " << score << " " << AlleleCountTypeToString(scoreToAlleleCountType(score)) << std::endl;

		std::lock_guard< std::mutex > l(m_counts_lock);
		auto iter = counts->find(samplePtr->getName());
		if (iter == counts->end())
		{
			std::vector< std::unordered_set< std::string > > sampleCounts((uint32_t)AlleleCountType::EndEnum);
			counts->emplace(samplePtr->getName(), sampleCounts);
			iter = counts->find(samplePtr->getName());
		}
		std::string bamID = bamAlignmentPtr->Name + std::to_string(bamAlignmentPtr->IsFirstMate());
		auto readCounts = iter->second;
		(*counts)[samplePtr->getName()][alleleCountType].emplace(bamID); // we are using readname so reads aren't counted more than once when we do the traceback and trackback through more than one reference node
	}

	std::unordered_set< std::string > Allele::getScoreCountFromAlleleCountType(const std::string& sampleName, AlleleCountType alleleCountType, bool forwardCount)
	{
		if (forwardCount)
		{
			auto iter = this->m_forward_counts.find(sampleName);
			if (iter != this->m_forward_counts.end())
			{
				return (iter->second)[(size_t)alleleCountType];
			}
		}
		else
		{
			auto iter = this->m_reverse_counts.find(sampleName);
			if (iter != this->m_reverse_counts.end())
			{
				return this->m_reverse_counts[sampleName][(size_t)alleleCountType];
			}
		}
		std::unordered_set< std::string > emptyValue;
		return emptyValue;
	}
}
