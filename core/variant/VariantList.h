#ifndef GRAPHITE_VARIANTLIST_H
#define GRAPHITE_VARIANTLIST_H

#include "IVariantList.h"
#include "core/region/Region.h"
#include "core/parser/VCFParser.hpp"

#include <vector>

namespace graphite
{
	class VariantList : public IVariantList
	{
	public:
		typedef std::shared_ptr< VariantList > SharedPtr;
		VariantList(const std::vector< IVariant::SharedPtr >& variantPtrs);
		~VariantList();

		bool getNextVariant(IVariant::SharedPtr& variantPtr) override;
		size_t getCount() override;
		void sort() override;
		void printToVCF(std::ostream& out, std::string& bamPath) override;
		void normalizeOverlappingVariants();
		void printHeader(std::ostream& out, std::string& bamPath);
		VariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr);
		void processOverlappingAlleles();

	protected:
		bool getNextCompoundVariant(IVariant::SharedPtr& variant);
		IVariant::SharedPtr buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< IVariant::SharedPtr >& variants);

		size_t m_current_index;
		std::vector< IVariant::SharedPtr > m_variant_ptrs;
		bool m_next_variant_init;
		IVariant::SharedPtr m_next_variant;
		VariantParser< const char* > m_vcf_parser;
	};
}

#endif //GRAPHITE_VARIANTLIST_H
