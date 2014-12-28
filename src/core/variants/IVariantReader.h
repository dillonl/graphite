#ifndef GWIZ_IVARIANT_READER_H
#define GWIZ_IVARIANT_READER_H

#include "core/variants/IVariant.h"
#include "core/Region.h"

#include "utils/NonCopyable.h"

#include <list>

namespace gwiz
{
	/*
	 * An iterface for variant readers. Provides a very sparse
	 * framework for a variant reader. This will allow flexibility
	 * for all file types and implementations.
	 */
	class IVariantReader : private noncopyable
	{
	public:
		typedef std::shared_ptr<IVariantReader> SharedPtr;
		IVariantReader() {}
		virtual ~IVariantReader() {}

		virtual std::list<IVariant::SharedPtr> GetAllVariants() = 0;
		virtual std::list<IVariant::SharedPtr> GetAllVariantsWithinChromosome(std::string& chrom) = 0;
		virtual std::list<IVariant::SharedPtr> GetVariantsWithinRegion(Region::SharedPtr region) = 0;
	};
} // end namespace gwiz

#endif //GWIZ_IVARIANT_READER_H
