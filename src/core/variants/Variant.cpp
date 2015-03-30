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
		for (auto variant : this->m_alt)
		{
			if (variant.size() > largest) { largest = variant.size(); }
		}
		return largest;
	}
}
