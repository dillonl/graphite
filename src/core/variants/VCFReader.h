#ifndef GWIZ_VCFREADER_H
#define GWIZ_VCFREADER_H

#include "IVariantReader.h"

#include <list>

namespace gwiz
{
	class VCFReader : public IVariantReader
	{
	public:
		typedef std::shared_ptr<VCFReader> SharedPtr;
		VCFReader(std::string& vcf_path);
		virtual ~VCFReader();

		virtual std::list<IVariant::SharedPtr> GetAllVariants();
		virtual std::list<IVariant::SharedPtr> GetAllVariantsWithinChromosome(std::string& chrom);
		virtual std::list<IVariant::SharedPtr> GetVariantsWithinRegion(Region::SharedPtr region);
	};
} // end namespace gwiz

#endif  //GWIZ_VCFREADER_H
