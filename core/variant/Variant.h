#ifndef GRAPHITE_VARIANT_H
#define GRAPHITE_VARIANT_H

#include <mutex>
#include <exception>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <tuple>
#include <stdlib.h>
#include <algorithm>

#include "IVariant.h"
#include "core/allele/Allele.h"
#include "core/reference/Reference.h"
#include "core/sample/Sample.h"
#include "core/util/Utility.h"

namespace graphite
{

	class IAlignment;
	class IAlignmentReader;

	class Variant : public IVariant
	{
	public:
		typedef std::shared_ptr< Variant > SharedPtr;
		typedef std::weak_ptr< Variant > WeakPtr;
		Variant(position pos, const std::string& chrom, const std::string& id, const std::string& quality, const std::string& filter, IAllele::SharedPtr refAllelePtr, std::vector< IAllele::SharedPtr > altAllelePtrs, uint32_t readLength=300);

		Variant();
		~Variant();

		static Variant::SharedPtr BuildVariant(const std::string& vcfLine, IReference::SharedPtr referencePtr, uint32_t readLength);
		void setSkip(bool skip) override;

		void setFilter(std::string filter) { this->m_filter = filter; }

		std::string getChrom() const override { return m_chrom; }
		position getPosition() override { return m_position; }
		std::string getQual() const { return m_qual; }
		std::string getFilter() const { return m_filter; }
		std::unordered_map< std::string, std::string > getInfoFields() const { return m_info_fields; }
		std::string getID() const { return m_id; }
		std::string getRef() { return std::string(this->m_ref_allele_ptr->getSequence()); }
		IAllele::SharedPtr getRefAllelePtr() override { return this->m_ref_allele_ptr; }
		std::vector< IAllele::SharedPtr > getAltAllelePtrs() override { return m_alt_allele_ptrs; }
		std::string getVariantLine(VCFHeader::SharedPtr headerPtr) override;
		void processOverlappingAlleles() override;

		uint32_t getAllelePrefixOverlapMaxCount(IAllele::SharedPtr allelePtr) override;
		uint32_t getAlleleSuffixOverlapMaxCount(IAllele::SharedPtr allelePtr) override;
		bool shouldSkip() override;

		std::vector< Region::SharedPtr > getRegions() override { return m_region_ptrs; }
		bool doesOverlap(IVariant::SharedPtr variantPtr, uint32_t readLength) override;
		uint32_t getReferenceSize() override;
		void addRegion(Region::SharedPtr regionPtr) override;
		uint32_t getVariantSize() override { return m_variant_size; }
		bool isStructuralVariant() { return m_is_sv; }

	protected:
		void setAlleleOverlapMaxCountIfGreaterThan(IAllele::SharedPtr allelePtr, std::unordered_map< IAllele::SharedPtr, uint32_t >& alleleOverlapCountMap, uint32_t overlapCount);
		void setRefAndAltAlleles(const std::string& ref, const std::vector< std::string >& alts);
		std::string getInfoFieldsString();
		std::string getSampleCounts(const std::string& sampleName);

		uint32_t m_max_prefix_match_length;
		uint32_t m_max_suffix_match_length;
		std::unordered_map< IAllele::SharedPtr, uint32_t > m_allele_prefix_max_overlap_map;
		std::unordered_map< IAllele::SharedPtr, uint32_t > m_allele_suffix_max_overlap_map;
		uint32_t m_position;
		std::string m_chrom;
		std::string m_id;
		std::string m_qual;
		std::string m_filter;
		std::string m_ref;
		std::vector< std::string > m_alt;
		IAllele::SharedPtr m_ref_allele_ptr;
		std::vector< IAllele::SharedPtr > m_alt_allele_ptrs;
		std::vector< IAllele::SharedPtr > m_all_allele_ptrs;
		std::unordered_map< std::string, std::string > m_info_fields;
		std::string m_line;
		bool m_skip;
		bool m_is_sv;
		uint32_t m_read_length;
		uint32_t m_reference_size;
		uint32_t m_variant_size;
		uint32_t m_overlap; // the amount the region overlaps the critical section of the variant (used for getting alignments from this region)

	private:
		void setUnorderedMapKeyValue(const std::string& rawString);
		void setAlleles(Reference::SharedPtr referencePtr, const std::string& vcfReferenceString, const std::string& alts, uint32_t readLength);
		void setAsDeletion(Reference::SharedPtr referencePtr, int svLength, uint32_t readLength);
		void setAsDuplication(Reference::SharedPtr referencePtr, int svLength, uint32_t readLength);
		void setAsInversion(Reference::SharedPtr referencePtr, int svLength, uint32_t readLength);
		void setAsInsertion(const std::string& ref, const std::string& alt, uint32_t readLength);
		void setAsStandardAlt(const std::string& ref, const std::string& alt, uint32_t readLength);

		std::string getFormatString();
		std::string getBlankSampleCounts();
		std::vector< Region::SharedPtr > m_region_ptrs;
	};

}

#endif //GRAPHITE_VARIANT_H
