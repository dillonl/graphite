#include "Variant.h"

namespace gwiz
{
	Variant::Variant()
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
			if (std::get< 1 >(allelePercentageTuple) > minAllowablePercentage)
			{
				vcfLines += std::get< 0 >(allelePercentageTuple);
				++printCount;
			}
			if (printCount >= maxAllelePrint) { break; }
		}
		return vcfLines;
	}
}
