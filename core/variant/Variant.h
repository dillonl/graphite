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
/*
#include "core/parser/VCFParser.hpp"
#include "core/parser/InfoFieldParser.hpp"
*/
#include "core/reference/Reference.h"
#include "core/alignment/Sample.hpp"
#include "core/util/Utility.h"

/*
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
*/

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

		inline static Variant::SharedPtr BuildVariant(const std::string& vcfLine, IReference::SharedPtr referencePtr, uint32_t maxAllowedAlleleSize=3000)
		{
			std::string fields;
			auto variantPtr = std::make_shared< Variant >();
			std::string ref;
			std::vector< std::string > alts;

			std::vector< std::string > vcfComponents;

			/* boost::split(vcfComponents, vcfLine, boost::is_any_of("\t")); */
			split(vcfLine, '\t', vcfComponents);

			if (vcfComponents.size() < 7)
			{
				std::cout << "vcf line is incorrectly formated" << std::endl;
			}

			variantPtr->m_chrom = vcfComponents[0];
			variantPtr->m_position = stoul(vcfComponents[1]);
			variantPtr->m_id = vcfComponents[2];
			ref = vcfComponents[3];
			if (vcfComponents[4].find(",") != std::string::npos)
			{
				/* boost::split(alts, vcfComponents[4], boost::is_any_of(",")); */
				split(vcfComponents[4], '\t', alts);
			}
			else
			{
				alts.push_back(vcfComponents[4]);
			}

			variantPtr->m_qual = vcfComponents[5];
			variantPtr->m_filter = vcfComponents[6];
			fields = vcfComponents[7];

			std::unordered_map< std::string, std::string > infoFields;
			setUnorderedMapKeyValue(fields, infoFields);
			/* static std::mutex mutLock; */
			/* std::lock_guard< std::mutex > lock(mutLock); */
			bool skipSequence = false;
			if (alts.size() == 1 &&  referencePtr != nullptr)
			{
				if (alts[0].compare("<DEL>") == 0)
				{
					int endPosition = stoi(infoFields["END"]);
					const char* reference = referencePtr->getSequence() + (variantPtr->m_position - referencePtr->getRegion()->getStartPosition());
					size_t alleleSize = endPosition - variantPtr->m_position;
					ref = std::string(reference, alleleSize);

					alts.clear();
					alts.emplace_back("");
				}
				else if (alts[0].compare("<DEP>") == 0)
				{
					try
					{
						int svLength = stoi(infoFields["SVLEN"]);

						const char* reference = referencePtr->getSequence() + (variantPtr->m_position - referencePtr->getRegion()->getStartPosition());
						size_t alleleSize = svLength / 2;
						size_t offset =  referencePtr->getRegion()->getStartPosition() - variantPtr->m_position;
						ref = std::string(reference, alleleSize);

						alts.clear();
						alts.emplace_back(ref + ref);
					}
					catch (int e)
					{
					}
					/* std::transform(refSequence.begin(), refSequence.end(),refSequence.begin(), ::toupper); */

				}
				else if (alts[0].compare("<INS>") == 0)
				{
					std::string altSequence = infoFields["SEQ"];
					std::transform(altSequence.begin(), altSequence.end(),altSequence.begin(), ::toupper);
					alts.clear();
					alts.emplace_back(altSequence);
					skipSequence = true;

				}
				else if (alts[0].compare("<INV>") == 0)
				{
					/* InfoFieldParser< std::string::iterator > infoParser; */
					/* std::map< std::string, std::string > infoFields; */
					//variantPtr->m_position -= 1;
					position offset = variantPtr->m_position - referencePtr->getRegion()->getStartPosition();
					auto svLen = abs(std::stoi(infoFields["SVLEN"])) - ref.size();
					/* std::cout << svLen << std::endl; */
					std::string altSequence = std::string(referencePtr->getSequence() + offset, svLen);
					std::reverse(altSequence.begin(), altSequence.end());
					/* std::cout << altSequence << std::endl; */
					alts.clear();
					alts.emplace_back(altSequence);
					/* position offset = variantPtr->m_position - referencePtr->getRegion()->getStartPosition(); */

					ref = std::string(referencePtr->getSequence() + offset, svLen);
				}
			}

			variantPtr->setRefAndAltAlleles(ref, alts);
			variantPtr->m_line = std::string(vcfLine.c_str());

			/* setUnorderedMapKeyValue(fields, variantPtr->m_info_fields); */
			variantPtr->setMaxAlleleSize();
			if (maxAllowedAlleleSize < variantPtr->m_max_allele_size && !skipSequence)
			{
				variantPtr->m_skip = true;
			}
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

		void setSkip(bool skip) override { m_skip = skip; }

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
		std::string getVariantLine(std::vector< std::shared_ptr< Sample > > samplePtrs) override;
		void processOverlappingAlleles() override;

		uint32_t getAllelePrefixOverlapMaxCount(IAllele::SharedPtr allelePtr) override;
		uint32_t getAlleleSuffixOverlapMaxCount(IAllele::SharedPtr allelePtr) override;
		void incrementUnmappedToMappedCount() override;
		void incrementMappedToUnmappedCount() override;
		void incrementRepositionedCount() override;
		bool shouldSkip() override;


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

		bool m_skip;

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
