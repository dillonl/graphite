#ifndef GWIZ_IGRAPH_H
#define GWIZ_IGRAPH_H

#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"
#include "core/variants/IVariant.h"
#include "core/graph/INode.h"

#include <boost/noncopyable.hpp>

#include <list>
#include <string>
#include <memory>

namespace gwiz
{
	class IGraph : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr<IGraph> SharedPtr;

	    IGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
		    m_reference_ptr(referencePtr), m_variant_list_ptr(variantListPtr), m_next_variant_init(false)
		{
		}
		virtual ~IGraph() {}

	protected:

		inline bool getNextCompoundVariant(Variant::SharedPtr& variant)
		{
			// the first time we call this function we need to get the first variant
			if (!m_next_variant_init)
			{
				m_next_variant_init = true;
				m_variant_list_ptr->getNextVariant(this->m_next_variant);
			}
			variant = this->m_next_variant; // set the variant to be returned
			if (this->m_next_variant == NULL) { return false; } // if the variant is NULL then return false, this indicates we are past the region of interest

			std::string referenceString = variant->getRef();
			std::vector< Variant::SharedPtr > variants;
			Variant::SharedPtr nextVariant;
			bool variantAdded = false;
			position startPosition = variant->getPosition();
			position compoundStartPosition = startPosition;
			position variantEndPosition = (startPosition + referenceString.size() - 1); // subtract 1 because we are counting starting with the position we are on

			// loop through variants until the variants stop overlapping.
			// As we loop through build a concatenated reference string
			// that represents the entire overlapped variants reference.
			// for those overlapped variants add them to a vector so
			// a "compound variant" can be generated.
			while (m_variant_list_ptr->getNextVariant(nextVariant))
			{
				if (variantEndPosition < nextVariant->getPosition())
				{
					break;
				}
				// this is a minor efficiency, even though this is a bit ugly
				// it is more efficient to add the variant after checking that this
				// is a compound variant because it is so rare
				if (!variantAdded)
				{
					variants.push_back(variant); // we will build a compound variant with all these variants
					variantAdded = true;
				}
				std::string nextReferenceString = nextVariant->getRef();
				position nextVariantEndPosition = (nextVariant->getPosition() + nextReferenceString.size()) - 1;
				// if the next variant has reference at a further position then add it to the end of the referenceString
				if (nextVariantEndPosition > variantEndPosition)
				{
					position referenceDelta = (nextVariantEndPosition - variantEndPosition);
					referenceString += nextReferenceString.substr(nextReferenceString.size() - referenceDelta);
					variantEndPosition = (compoundStartPosition + referenceString.size() - 1); // subtract 1 because we are counting starting with the position we are on
				}
				variants.push_back(nextVariant); // we will build a compound variant with all these variants
				if (nextVariant->getPosition() < startPosition)
				{
					startPosition = nextVariant->getPosition();
				}
			}
			if (!variants.empty())
			{
				variant = buildCompoundVariant(startPosition, referenceString, variants);
			}
			this->m_next_variant = nextVariant; // set the next variant
			return true;
		}

		Variant::SharedPtr buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< Variant::SharedPtr >& variants)
		{
			/*
			if (true)
			{
				return variants[0];
			}
			*/
			position referenceEndPosition = referenceString.size() + startPosition;
			std::string chrom = variants[0]->getChrom();
			position pos = variants[0]->getPosition();
			std::string id = ".";
			std::string line = chrom + "\t" + std::to_string(pos) + "\t" + id + "\t" + referenceString + "\t";
			std::map< std::string, std::string > altMap; // maps new alt with the original vcf line
			// loop over all the variants
			for (auto variantPtr : variants)
			{
				// loop over all the alts in the variants
				for (auto altString : variantPtr->getAlt())
				{
					// basically we are replacing the variant's reference with the alt within the aggrigated reference (referenceString)
					// it's complicated to explain in words but if you follow the code it isn't too bad
					std::string variantString = referenceString;
					variantString.erase(variantPtr->getPosition() - startPosition, variantPtr->getRef().size());
					variantString.insert(variantPtr->getPosition() - startPosition, altString);
					line += variantString + ",";
					altMap[variantString] = variantPtr->getVCFLineFromAlternate(altString);
				}
			}
			/* line.replace(line.size() - 1, 1, "\t\n"); // replace the past comma with a tab */
			line.replace(line.size() - 1, 1, "\t"); // replace the past comma with a tab
			std::string quality = "0";
			std::string filter = "PASS";
			std::string info = "C=TRUE";
			line += quality + "\t" + filter + "\t" + info + "\t";
			auto variant = Variant::BuildVariant(line.c_str(), this->m_vcf_parser);
			for (auto const& altTuple : altMap) // set the vcf_lines back to the original vcf line
			{
				variant->setVCFLineFromAlternate(altTuple.first, altTuple.second.c_str(), altTuple.second.size());
			}
			return variant;
		}

		virtual void constructGraph() = 0;

		IReference::SharedPtr m_reference_ptr;
		IVariantList::SharedPtr m_variant_list_ptr;

		// This is used in conjunction with getCompoundNode.
		// This stores the next variant so we don't have to
		// get the variant twice
		Variant::SharedPtr m_next_variant;
		VariantParser< const char* > m_vcf_parser;
	private:
		bool m_next_variant_init;


	};
} // end namespace gwiz

#endif // GWIZ_IGRAPH_H
