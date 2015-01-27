#ifndef GWIZ_IVARIANT_READER_H
#define GWIZ_IVARIANT_READER_H

#include "core/variants/IVariant.h"
#include "core/variants/Variant.h"

#include "core/region/Region.h"

#include "core/utils/NonCopyable.h"

#include <list>

namespace gwiz
{
	/*
	 * An iterface for variant readers. Provides a very sparse
	 * framework for a variant reader. This will allow flexibility
	 * for all file types and implementations.
	 */
	class IVariantList : private noncopyable
	{
	public:
		typedef std::shared_ptr< IVariantList > SharedPtr;
		IVariantList() {}
		virtual ~IVariantList() {}

		virtual inline bool getNextVariant(Variant::SharedPtr& variant) = 0;
	};
} // end namespace gwiz

#endif //GWIZ_IVARIANT_READER_H
