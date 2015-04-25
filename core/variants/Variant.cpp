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
		alleleCountString += std::to_string(this->m_allele_count[getRef()]) + ",";
		alleleCountString += std::to_string(this->m_allele_reverse_strand_count[getRef()]);
		for (auto& alt : getAlt())
		{
			alleleCountString += ",";
			alleleCountString += std::to_string(this->m_allele_count[alt]) + ",";
			alleleCountString += std::to_string(this->m_allele_reverse_strand_count[alt]);
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
			auto arc = this->m_allele_reverse_strand_count.find(alt);
			if ((ac != this->m_allele_count.end() && ac->second > 0) ||
				(arc != this->m_allele_reverse_strand_count.end() && arc->second > 0))
			{
				return true;
			}
		}
		return false;
	}

	void Variant::printVariant(std::ostream& out)
	{
		std::string passFail = (this->m_pass) ? "PASS" : "FAIL";
		out << this->m_chrom << "\t" << getPosition() << "\t.\t" << alleleString() << "\tNA\t" << passFail << "\tDP=" << this->m_total_allele_count << ";DP4" << getAlleleCountString() << std::endl;
	}
}
