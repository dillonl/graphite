#include "Variant.h"

namespace gwiz
{
	Variant::Variant() :
		m_total_allele_count(0)
	{
	}

	Variant::~Variant()
	{
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
		alleleCountString += "R+:" + std::to_string(this->m_allele_count[getRef()]) + " ";
		alleleCountString += "R-:" + std::to_string(this->m_allele_reverse_strand_count[getRef()]) + " ";
		size_t count = 1;
		for (auto& alt : getAlt())
		{
			alleleCountString += "A" + std::to_string(count) + "+:" + std::to_string(this->m_allele_count[alt]) + " ";
			alleleCountString += "A" + std::to_string(count) + "-:" + std::to_string(this->m_allele_reverse_strand_count[alt]) + " ";
			++count;
		}
		alleleCountString.pop_back(); // removes the last space off the end
		return alleleCountString;
	}

	std::string Variant::toString()
	{
		std::string vcfLines = "";
		std::vector< std::tuple< std::string, uint32_t > > allelePercentages;
		generateAlleleCoveragePercentages(allelePercentages);
		uint32_t printCount = 0;
		uint32_t maxAllelePrint = 2; // the max number of vcf lines to print; 2 for a diploid
		uint32_t minAllowablePercentage = 30; // the minimum percentage acceptable before printing
		for (auto& allelePercentageTuple : allelePercentages)
		{
			// std::cout << std::get< 0 >(allelePercentageTuple) << " " << std::get< 1 >(allelePercentageTuple) << std::endl;
			if (std::get< 1 >(allelePercentageTuple) > minAllowablePercentage)
			{
				vcfLines += std::to_string(getPosition()) + " " + getAlleleCountString() + "\n";
				// vcfLines += getVCFLineFromAlternate(std::get< 0 >(allelePercentageTuple));
				++printCount;
			}
			if (printCount >= maxAllelePrint) { break; } // only need to print 2 (diploid)
		}
		return vcfLines;
	}
}
