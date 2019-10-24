#include "Allele.h"

#include "core/graph/Node.h"

namespace graphite
{
	Allele::Allele(const std::string& sequence) :
		m_sequence(sequence)
	{
	}

	Allele::~Allele()
	{
	}

	void Allele::incrementScoreCount(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, int score)
	{
		bool isForwardStrand = !bamAlignmentPtr->IsReverseStrand();
		// static std::mutex slock;
		// std::lock_guard< std::mutex > lg(slock);
		size_t alleleCountType = (size_t)scoreToAlleleCountType(score);
		if (score < 0)
		{
			alleleCountType = (size_t)AlleleCountType::Ambiguous;
		}
		// std::cout << "ambiguous type: " << AlleleCountTypeToString((AlleleCountType)alleleCountType) << std::endl;
		std::unordered_map< std::string, std::vector< std::unordered_set< std::string > > >* counts = (isForwardStrand) ? &this->m_forward_counts : &this->m_reverse_counts;
		for (auto allelePtr : this->m_paired_allele_ptrs)
		{
			allelePtr->incrementScoreCount(bamAlignmentPtr, samplePtr, score);
		}

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

	void Allele::pairAllele(Allele::SharedPtr allelePtr)
	{
		this->m_paired_allele_ptrs.emplace(allelePtr);
	}

	void Allele::addSemanticLoci(position pos, const std::string& refSequence, const std::string& altSequence)
	{
		auto iter = this->m_semantic_locations.find(pos);
		if (iter == this->m_semantic_locations.end())
		{
			std::unordered_set< std::string > sequences = {refSequence + ":" + altSequence};
			this->m_semantic_locations.emplace(pos, sequences);
		}
		else
		{
			iter->second.emplace(refSequence + ":" + altSequence);
		}
	}

	void Allele::registerSupportingReadInformation(SupportingReadInfo::SharedPtr supportingReadInfo)
	{
		std::lock_guard< std::mutex > l(this->m_supporting_read_info_mutex);
		this->m_supporting_read_info_list.emplace_back(supportingReadInfo);
	}

	std::vector< SupportingReadInfo::SharedPtr > Allele::getSupportingReadInfoPtrs()
	{
		std::lock_guard< std::mutex > l(this->m_supporting_read_info_mutex);
		std::vector< SupportingReadInfo::SharedPtr > supportingReadInfoList(this->m_supporting_read_info_list);
		return supportingReadInfoList;
	}
}
