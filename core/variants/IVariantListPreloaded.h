#ifndef GWIZ_IVARIANTLISTPRELOADED_H
#define GWIZ_IVARIANTLISTPRELOADED_H

#include "IVariantList.h"

namespace gwiz
{
	class IVariantListPreloaded : public IVariantList
	{
	public:
		typedef std::shared_ptr< IVariantListPreloaded > SharedPtr;
		IVariantListPreloaded() {}
		virtual ~IVariantListPreloaded() {}

		virtual bool getPreviousVariant(Variant::SharedPtr& variant) = 0;
		virtual bool getVariant(Variant::SharedPtr& variant, const uint32_t index) = 0;
	};
}

#endif //GWIZ_IVARIANTLISTPRELOADED_H
