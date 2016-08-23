#ifndef GRAPHITE_VARIANTLIST_H
#define GRAPHITE_VARIANTLIST_H

#include "IVariantList.h"
#include "core/region/Region.h"
#include "core/parser/VCFParser.hpp"
#include "core/reference/IReference.h"
#include "VCFHeader.h"

#include <zlib.h>

#include <vector>

namespace graphite
{
	class VariantList : public IVariantList
	{
	public:
		typedef std::shared_ptr< VariantList > SharedPtr;
		VariantList(const std::vector< IVariant::SharedPtr >& variantPtrs, IReference::SharedPtr referencePtr);
		~VariantList();

		bool getNextVariant(IVariant::SharedPtr& variantPtr) override;
		size_t getCount() override;
		void sort() override;
		void printToVCF(IHeader::SharedPtr header, bool printHeader, std::ostream& out) override;
		void normalizeOverlappingVariants();
		void printHeader(std::ostream& out, std::string& bamPath);
		void printToCompressedVCF(IHeader::SharedPtr headerPtr, bool printHeader, int out);
		VariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr);
		void processOverlappingAlleles() override;

	protected:
		bool getNextCompoundVariant(IVariant::SharedPtr& variant);
		IVariant::SharedPtr buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< IVariant::SharedPtr >& variants);

		size_t m_current_index;
		std::vector< IVariant::SharedPtr > m_variant_ptrs;
		bool m_next_variant_init;
		IVariant::SharedPtr m_next_variant;
		VariantParser< const char* > m_vcf_parser;
		IReference::SharedPtr m_reference_ptr;
	};
}

#endif //GRAPHITE_VARIANTLIST_H
