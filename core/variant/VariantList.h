#ifndef GRAPHITE_VARIANTLIST_H
#define GRAPHITE_VARIANTLIST_H

#include "core/region/Region.h"
#include "core/reference/IReference.h"
#include "VCFHeader.h"
#include "IVariant.h"
#include "core/file/IFileWriter.h"
#include "core/util/Noncopyable.hpp"

#include <zlib.h>

#include <vector>

namespace graphite
{
	class VariantList : private Noncopyable
	{
	public:
		typedef std::shared_ptr< VariantList > SharedPtr;
		VariantList(const std::vector< IVariant::SharedPtr >& variantPtrs, IReference::SharedPtr referencePtr);
		~VariantList();

		bool getNextVariant(IVariant::SharedPtr& variantPtr);
		bool peekNextVariant(IVariant::SharedPtr& variantPtr);
		size_t getCount();
		void sort();
		void normalizeOverlappingVariants();
		void printHeader(std::ostream& out, std::string& bamPath);
		void printToCompressedVCF(VCFHeader::SharedPtr headerPtr, bool printHeader, int out);
		VariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr);
		void processOverlappingAlleles();
		void writeVariantList(IFileWriter::SharedPtr fileWriter, VCFHeader::SharedPtr headerPtr, bool printHeader);
		std::vector< IVariant::SharedPtr > getAllVariantPtrs() { return m_variant_ptrs; }

	protected:
		bool getNextCompoundVariant(IVariant::SharedPtr& variant);
		IVariant::SharedPtr buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< IVariant::SharedPtr >& variants);

		size_t m_current_index;
		std::vector< IVariant::SharedPtr > m_variant_ptrs;
		bool m_next_variant_init;
		IVariant::SharedPtr m_next_variant;
		IReference::SharedPtr m_reference_ptr;
	};
}

#endif //GRAPHITE_VARIANTLIST_H
