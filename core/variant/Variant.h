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

#include "IVariant.h"
#include "core/allele/Allele.h"
#include "core/parser/VCFParser.hpp"
#include "core/parser/InfoFieldParser.hpp"
#include "core/reference/Reference.h"
#include "core/alignment/Sample.hpp"

#include <boost/algorithm/string.hpp>

namespace graphite
{

	class IAlignment;
	class IAlignmentReader;
	class Variant : public IVariant
	{
	public:
		typedef std::shared_ptr< Variant > SharedPtr;
		typedef std::weak_ptr< Variant > WeakPtr;
		Variant(position pos, const std::string& chrom, const std::string& id, const std::string& quality, const std::string& filter, IAllele::SharedPtr refAllelePtr, std::vector< IAllele::SharedPtr > altAllelePtrs);

		Variant();
		~Variant();

		inline static Variant::SharedPtr BuildVariant(const std::string& vcfLine, VariantParser< const char* >& parser, IReference::SharedPtr referencePtr)
		{
			std::string fields;
			auto variantPtr = std::make_shared< Variant >();
			std::string ref;
			std::vector< std::string > alts;
			if (!boost::spirit::qi::parse(vcfLine.c_str(), (vcfLine.c_str() + vcfLine.size()), parser, variantPtr->m_chrom, variantPtr->m_position, variantPtr->m_id, ref, alts, variantPtr->m_qual, variantPtr->m_filter, fields))
			{
				throw "An invalid line in the VCF caused an exception. Please correct the input and try again";
			}

			if (alts.size() == 1 && alts[0].compare("<DEL>") == 0 && referencePtr != nullptr)
			{
				InfoFieldParser< std::string::iterator > infoParser;
				std::map< std::string, std::string > infoFields;
				if (boost::spirit::qi::parse(fields.begin(), fields.end(), infoParser, infoFields))
				{
					variantPtr->m_position -= 1;
					position offset = variantPtr->m_position - referencePtr->getRegion()->getStartPosition();
					ref = std::string(referencePtr->getSequence() + offset, abs(std::stoi(infoFields["SVLEN"])) + 1);
					alts[0] = std::string(referencePtr->getSequence() + offset, 1);
				}
			}

			variantPtr->setRefAndAltAlleles(ref, alts);
            variantPtr->m_line = std::string(vcfLine.c_str(), vcfLine.size() - 1);

			/* setUnorderedMapKeyValue(fields, variantPtr->m_info_fields); */
			variantPtr->setMaxAlleleSize();
			return variantPtr;
		}

		inline static void setUnorderedMapKeyValue(const std::string& rawString, std::unordered_map< std::string, std::string >& keyValue)
		{
			const char* raw = rawString.c_str();
			std::string key;
			std::string value;
			char prevToken = ';';
			for (size_t i = 0; i < rawString.size(); ++i)
			{
				if (raw[i] == '=')
				{
					prevToken = '=';
				}
				else if (raw[i] == ';')
				{
					keyValue[key] = value;
					value = "";
					key = "";
					prevToken = ';';
				}
				else if (prevToken == '=')
				{
					value += raw[i];
				}
				else if (prevToken == ';')
				{
					key += raw[i];
				}
			}
			keyValue[key] = value;
		}

		inline bool processSV(IReference::SharedPtr referencePtr)
		{
			auto cn0AltIter = std::find(this->m_alt.begin(), this->m_alt.end(), "<CN0>");
			auto endInfoIter = this->m_info_fields.find("END");
			std::string endResult = "none";
			/*
			for (auto iter : this->m_info_fields)
			{
				std::cout << iter.first << " " << iter.second << std::endl;
			}
			exit(0);
			*/
			/* std::string endExists = (endInfoIter == this->m_info_fields.end()) ? "false" : "true"; */
			/* std::cout << this->m_position << " " << endResult << " "  << endExists; */
			if (cn0AltIter != this->m_alt.end() && endInfoIter != this->m_info_fields.end())
			{
				position endPosition = atoi(endInfoIter->second.c_str());
				if (endPosition > referencePtr->getRegion()->getEndPosition()) { return false; }
				/* std::string endPosition = this->m_info_fields["END"]; */
				this->m_alt.clear();
				this->m_alt.push_back(this->getRef());
				const char* reference = referencePtr->getSequence() + (this->m_position - referencePtr->getRegion()->getStartPosition());
				size_t alleleSize = endPosition - this->m_position;
				size_t offset =  referencePtr->getRegion()->getStartPosition() - this->m_position;
				/* std::cout << m_position << std::endl; */
				/*
				std::cout << "-------" << std::endl;
				std::cout << referenceAllele << std::endl;
				std::cout << "-------" << std::endl;
				*/
				this->m_ref = std::string(reference, alleleSize);
			}
			else
			{
				for (auto& alt : this->m_alt)
				{
					/* if (alt.c_str()[0] == '<') { return false; } */
				}
			}
			return true;
		}

		void setFilter(std::string filter) { this->m_filter = filter; }
		/* std::string getAlleleCountString(); */
		std::string alleleString();
		/* void addPotentialAlignment(const std::shared_ptr< IAlignment > alignmentPtr); */

		std::string getChrom() const override { return m_chrom; }
		position getPosition() override { return m_position; }
		std::string getQual() const { return m_qual; }
		std::string getFilter() const { return m_filter; }
		std::unordered_map< std::string, std::string > getInfoFields() const { return m_info_fields; }
		std::string getID() const { return m_id; }
		std::string getRef() { return std::string(this->m_ref_allele_ptr->getSequence()); }
		size_t getMaxAlleleSize() override { return this->m_max_allele_size; }
		IAllele::SharedPtr getRefAllelePtr() override { return this->m_ref_allele_ptr; }
		std::vector< IAllele::SharedPtr > getAltAllelePtrs() override { return m_alt_allele_ptrs; }
		void printVariant(std::ostream& out, std::vector< std::shared_ptr< Sample > > samplePtrs) override;
		void processOverlappingAlleles() override;

		uint32_t getAllelePrefixOverlapMaxCount(IAllele::SharedPtr allelePtr) override;
		uint32_t getAlleleSuffixOverlapMaxCount(IAllele::SharedPtr allelePtr) override;
		void incrementUnmappedToMappedCount() override;
		void incrementMappedToUnmappedCount() override;
		void incrementRepositionedCount() override;

	protected:
		void setAlleleOverlapMaxCountIfGreaterThan(IAllele::SharedPtr allelePtr, std::unordered_map< IAllele::SharedPtr, uint32_t >& alleleOverlapCountMap, uint32_t overlapCount);
		void setRefAndAltAlleles(const std::string& ref, const std::vector< std::string >& alts);
		std::string getGenotype();
		std::string getInfoFieldsString();
		/* uint32_t getTotalAlleleCount(); */
		std::string getSampleCounts(std::vector< Sample::SharedPtr > samplePtrs);

		size_t m_max_allele_size;
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
		std::string m_line;
		IAllele::SharedPtr m_ref_allele_ptr;
		std::vector< IAllele::SharedPtr > m_alt_allele_ptrs;
		std::vector< IAllele::SharedPtr > m_all_allele_ptrs;
		std::unordered_map< std::string, std::string > m_info_fields;
		std::atomic< uint32_t > m_unmapped_to_mapped_count;
		std::atomic< uint32_t > m_mapped_to_unmapped_count;
		std::atomic< uint32_t > m_repositioned_count;

	private:
		void setMaxAlleleSize();

		// a helper class for the getSampleCounts method
		class VariantSampleContainer
		{
		public:
		    VariantSampleContainer() :
			    m_total_count(0),
				m_forward_count(0),
				m_reverse_count(0)
			{
			}

			~VariantSampleContainer()
			{
			}

			uint32_t getTotalCount() { return this->m_total_count; }
			uint32_t getForwardCount() { return this->m_forward_count; }
			uint32_t getReverseCount() { return this->m_reverse_count; }

			void addToTotalCount(uint32_t value) { this->m_total_count += value; }
			void addToForwardCount(uint32_t value) { this->m_forward_count += value; }
			void addToReverseCount(uint32_t value) { this->m_reverse_count += value; }
		private:
			uint32_t m_total_count;
			uint32_t m_forward_count;
			uint32_t m_reverse_count;
		};
	};

}

#endif //GRAPHITE_VARIANT_H
