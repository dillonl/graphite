#ifndef GWIZ_VARIANTLISTVCFPRELOADED_H
#define GWIZ_VARIANTLISTVCFPRELOADED_H

#include "VariantList.h"

namespace gwiz
{
	class VariantListVCFPreloaded : public VariantList
	{
	public:
		typedef std::shared_ptr< VariantListVCFPreloaded > SharedPtr;
		VariantListVCFPreloaded(const std::string& path, Region::SharedPtr region);
		VariantListVCFPreloaded(const std::string& path);
		~VariantListVCFPreloaded();

		/* bool getNextVariant(Variant::SharedPtr& variant) override; */
		/* bool getPreviousVariant(Variant::SharedPtr& variant) override; */
		/* bool getVariant(Variant::SharedPtr& variant, const uint32_t index) override; */

		void loadVariantsFromFile();

	protected:
		std::vector< Variant::SharedPtr > m_variants;
		std::string m_path;
		Region::SharedPtr m_region_ptr;
		bool m_processed;
	};
}

#endif //GWIZ_VARIANTLISTVCFPRELOADED_H
