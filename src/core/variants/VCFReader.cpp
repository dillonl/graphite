#include "VCFReader.h"

namespace gwiz
{
	VCFReader::VCFReader(std::string& vcf_path)
	{
	}

	VCFReader::~VCFReader()
	{
	}

	std::list<IVariant::SharedPtr> VCFReader::GetAllVariants()
	{
		std::list<IVariant::SharedPtr> variants;
		return variants;
	}

	std::list<IVariant::SharedPtr> VCFReader::GetAllVariantsWithinChromosome(std::string& chrom)
	{
		std::list<IVariant::SharedPtr> variants;
		return variants;
	}

	std::list<IVariant::SharedPtr> VCFReader::GetVariantsWithinRegion(Region::SharedPtr region)
	{
		std::list<IVariant::SharedPtr> variants;
		return variants;
	}
} // end namespace gwiz
