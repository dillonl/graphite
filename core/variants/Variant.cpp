#include "Variant.h"
#include "core/alignments/IAlignment.h"

namespace gwiz
{
	Variant::Variant() :
		m_total_allele_count_low_quality(0),
		m_total_allele_count(0)
	{
	}

	Variant::~Variant()
	{
	}

	void Variant::initializeAlleleCounters()
	{
		m_allele_count[getRef()] = std::make_tuple(0, 0);
		for (const auto& alt : getAlt())
		{
			m_allele_count[alt] = std::make_tuple(0, 0);
		}
	}

	size_t Variant::getSmallestAlleleSize()
	{
		size_t smallest = this->m_ref[0].size();
		for (auto variant : this->m_alt)
		{
			if (variant.size() < smallest) { smallest = variant.size(); }
		}
		return smallest;
	}

	size_t Variant::getLargestAlleleSize()
	{
		size_t largest = this->m_ref[0].size();
		for (auto& variant : this->m_alt)
		{
			if (variant.size() > largest) { largest = variant.size(); }
		}
		return largest;
	}

	std::string Variant::getAlleleCountString()
	{
		std::string alleleCountString = "";
		auto refAlleleTuple = this->m_allele_count[getRef()];
		alleleCountString += std::to_string(std::get< 0 >(refAlleleTuple)) + ",";
		alleleCountString += std::to_string(std::get< 1 >(refAlleleTuple));
		for (auto& alt : getAlt())
		{
			auto altAlleleTuple = this->m_allele_count[alt];
			alleleCountString += ",";
			alleleCountString += std::to_string(std::get< 0 >(altAlleleTuple)) + ",";
			alleleCountString += std::to_string(std::get< 1 >(altAlleleTuple));
		}
		return alleleCountString;
	}

	std::string Variant::alleleString()
	{
		std::string alleleString = "";
		alleleString += getRef() + "\t";
		for (auto& alt : getAlt())
		{
			alleleString += alt + ",";
		}
		alleleString.pop_back(); // removes the last space off the end
		return alleleString;
	}

	bool Variant::hasAlts()
	{
		for (auto& alt : getAlt())
		{
			auto ac = this->m_allele_count.find(alt);
			if ((ac != this->m_allele_count.end() && std::get< 0 >(ac->second) > 0) ||
				(ac != this->m_allele_count.end() && std::get< 1 >(ac->second) > 0))
			{
				return true;
			}
		}
		return false;
	}

	void Variant::addPotentialAlignment(const IAlignment::SharedPtr alignmentPtr, const std::string& allele)
	{
		this->m_potential_alignments.emplace(std::make_pair(alignmentPtr, allele));
	}

	void Variant::calculateAlleleCounts()
	{
		for (const auto& alignmentPairIter : this->m_potential_alignments)
		{
			int32_t mappingScore = alignmentPairIter.first->getVariantMappingScore(getVariantID());
			if (0 <= mappingScore)
			{
				if (alignmentPairIter.first->isMapped())
				{
					increaseCount(alignmentPairIter.second, alignmentPairIter.first);
				}
				else
				{
					incrementLowQualityCount(alignmentPairIter.first);
				}
			}
		}
	}

	void Variant::printVariant(std::ostream& out)
	{
		calculateAlleleCounts();
		uint32_t totalCount = this->m_total_allele_count + this->m_total_allele_count_low_quality;
		out << this->m_chrom << "\t" << getPosition() << "\t.\t" << alleleString() << "\t0\t.\tDP=" << this->m_total_allele_count << ";DP4=" << getAlleleCountString() << ";TC=" << totalCount << std::endl;
	}

}
