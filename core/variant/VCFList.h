#ifndef GWIZ_VCFLIST_H
#define GWIZ_VCFLIST_H

#include "IVariantList.h"
#include "Variant.h"
#include "core/parser/VCFParser.hpp"
#include "core/region/Region.h"

#include <vector>

namespace gwiz
{
	/*
	 * The VCFList is a IVariantList that has been read from
	 * a file. This has very limited functionality due to the
	 * list being read from the file rather than a container.
	 */
	class VCFList : public IVariantList
	{
	public:
		typedef std::shared_ptr< VCFList > SharedPtr;
		VCFList(const std::string& filePath, Region::SharedPtr regionPtr);
		~VCFList();

		bool getNextVariant(IVariant::SharedPtr& variantPtr) override;
		void printToVCF(std::ostream& out) override;
		bool hasVariants() override;
		std::vector< IVariant::SharedPtr > load();

		//things to be removed possibly
		/* size_t getCount() override; */
		/* void sort() override; */
		/* void addVariants(IVariantList::SharedPtr variantListPtr) override; */

	private:

	};
}

#endif //GWIZ_VCFLIST_H
