#ifndef GWIZ_VCFVARIANT_H
#define GWIZ_VCFVARIANT_H

#include "Variant.h"
#include "VCFFileReader.h"

namespace gwiz
{
	class VCFVariant : public Variant
	{
	public:
		~VCFVariant() {}

		inline static VCFVariant::SharedPtr BuildVariant(const std::string& vcfLine, VCFFileReader::SharedPtr vcfReaderPtr, VariantParser< const char* >& parser)
		{
			auto variantPtr = std::make_shared< VCFVariant >();

			return variantPtr;
		}

		VCFVariant() {}

	protected:
	};
}

#endif //GWIZ_VCFVARIANT_H
