#ifndef GRAPHITE_VCFVARIANT_H
#define GRAPHITE_VCFVARIANT_H

#include "Variant.h"
#include "VCFFileReader.h"

namespace graphite
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

#endif //GRAPHITE_VCFVARIANT_H
