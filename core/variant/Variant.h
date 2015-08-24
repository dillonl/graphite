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

#include "IVariant.h"
#include "core/allele/Allele.h"
#include "core/parser/VCFParser.hpp"
#include "core/reference/Reference.h"


namespace graphite
{

	class IAlignment;
	class Variant : public IVariant
	{
	public:
		typedef std::shared_ptr< Variant > SharedPtr;
		Variant(position pos, const std::string& chrom, const std::string& id, const std::string& quality, const std::string& filter, IAllele::SharedPtr refAllelePtr, std::vector< IAllele::SharedPtr > altAllelePtrs);

		Variant();
		~Variant();

		inline static Variant::SharedPtr BuildVariant(const std::string& vcfLine, VariantParser< const char* >& parser)
		{
			std::string fields;
			auto variantPtr = std::make_shared< Variant >();
			std::string ref;
			std::vector< std::string > alts;
			if (!boost::spirit::qi::parse(vcfLine.c_str(), (vcfLine.c_str() + vcfLine.size()), parser, variantPtr->m_chrom, variantPtr->m_position, variantPtr->m_id, ref, alts, variantPtr->m_qual, variantPtr->m_filter, fields))
			{
				throw "An invalid line in the VCF caused an exception. Please correct the input and try again";
			}

			variantPtr->setRefAllele(ref);
			variantPtr->setAltAlleles(alts);

			variantPtr->m_all_allele_ptrs.reserve(variantPtr->m_alt_allele_ptrs.size() + 1);
			variantPtr->m_all_allele_ptrs.emplace_back(variantPtr->m_ref_allele_ptr);
			variantPtr->m_all_allele_ptrs.insert(variantPtr->m_all_allele_ptrs.end(), variantPtr->m_alt_allele_ptrs.begin(), variantPtr->m_alt_allele_ptrs.end());

			setUnorderedMapKeyValue(fields, variantPtr->m_info_fields);
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

		void setFilter(std::string filter)
		{
			this->m_filter = filter;
		}

		std::string getAlleleCountString();
		std::string alleleString();
		/* void addPotentialAlignment(const std::shared_ptr< IAlignment > alignmentPtr); */

		std::string getChrom() const { return m_chrom; }
		position getPosition() override { return m_position; }
		std::string getQual() const { return m_qual; }
		std::string getFilter() const { return m_filter; }
		std::unordered_map< std::string, std::string > getInfoFields() const { return m_info_fields; }
		std::string getID() const { return m_id; }
		std::string getRef() { return std::string(this->m_ref_allele_ptr->getSequence()); }
		IAllele::SharedPtr getRefAllelePtr() override { return this->m_ref_allele_ptr; }
		std::vector< IAllele::SharedPtr > getAltAllelePtrs() override { return m_alt_allele_ptrs; }
		void printVariant(std::ostream& out) override;
		void processOverlappingAlleles() override;

	protected:
		static std::vector< std::string > removeDuplicateAltAlleles(const std::vector< std::string >& alts)
		{
			std::vector< std::string > tmpAlts;
			std::unordered_map< std::string, bool > seqMap;
			for (auto alt : alts)
			{
				if (seqMap.find(alt) == seqMap.end())
				{
					tmpAlts.emplace_back(alt);
					seqMap[alt] = true;
				}
			}
			return tmpAlts;
		}

		void setRefAllele(const std::string& ref)
		{
			this->m_ref_allele_ptr = std::make_shared< Allele >(SequenceManager::Instance()->getSequence(ref.c_str()));
		}

		void setAltAlleles(const std::vector< std::string >& alts)
		{
			this->m_alt_allele_ptrs.clear();
			for (const auto& alt : alts)
			{
				auto sequencePtr = SequenceManager::Instance()->getSequence(alt.c_str());
				auto altAllelePtr = std::make_shared< Allele >(sequencePtr);
				this->m_alt_allele_ptrs.emplace_back(altAllelePtr);
			}
		}

		std::string getGenotype();
		uint32_t getTotalAlleleCount();

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
	};

}

#endif //GRAPHITE_VARIANT_H
