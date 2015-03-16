#ifndef GWIZ_VARIANTLISTVCFPRELOADED_H
#define GWIZ_VARIANTLISTVCFPRELOADED_H

#include "IVariantListPreloaded.h"

namespace gwiz
{
	class VariantListVCFPreloaded : public IVariantListPreloaded //, public VCFFileReader
	{
	public:
		typedef std::shared_ptr< VariantListVCFPreloaded > SharedPtr;
		VariantListVCFPreloaded(const std::string& path, Region::SharedPtr region);
		VariantListVCFPreloaded(const std::string& path);
		~VariantListVCFPreloaded();

		bool getNextVariant(Variant::SharedPtr& variant) override;
		bool getPreviousVariant(Variant::SharedPtr& variant) override;
		bool getVariant(Variant::SharedPtr& variant, const uint32_t index) override;

		void process();

	protected:
		std::vector< Variant::SharedPtr > m_variants;
		bool m_processed;
	};
}

#endif //GWIZ_VARIANTLISTVCFPRELOADED_H
