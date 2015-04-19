#ifndef GWIZ_IVARIANTLIST_H
#define GWIZ_IVARIANTLIST_H

#include <boost/noncopyable.hpp>

#include "core/variants/IVariant.h"
#include "core/variants/Variant.h"

#include "core/region/Region.h"

#include <vector>
#include <algorithm>

namespace gwiz
{
	/*
	 * An iterface for variant readers. Provides a very sparse
	 * framework for a variant reader. This will allow flexibility
	 * for all file types and implementations.
	 */
	class IVariantList : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IVariantList > SharedPtr;
		IVariantList() {}
		virtual ~IVariantList() {}

		virtual bool getNextVariant(Variant::SharedPtr& variant) = 0;
		virtual void rewind() = 0;
		virtual IVariantList::SharedPtr getVariantsInRegion(Region::SharedPtr region) = 0;
		virtual void addVariants(IVariantList::SharedPtr variantsListPtr) {}
		virtual void addVariant(Variant::SharedPtr variantPtr) {}
		virtual void sortVariants() {}
		virtual size_t getCount() = 0;

		virtual void printToVCF(std::ostream& out) = 0;
	protected:

	};
} // end namespace gwiz

#endif //GWIZ_IVARIANTLIST_H
