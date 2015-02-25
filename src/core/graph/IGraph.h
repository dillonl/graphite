#ifndef GWIZ_IGRAPH_H
#define GWIZ_IGRAPH_H

#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"
#include "core/variants/IVariant.h"
#include "core/graph/INode.h"
#include "../plugins/vg/graph/ReferenceNode.h"

#include <boost/noncopyable.hpp>

#include <list>
#include <string>
#include <memory>

namespace gwiz
{
	template <typename T>
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
		IGraph() {} // used in tests

		virtual T addVariantNode(INode::SharedPtr variantNodePtr);
		virtual T addReferenceNodeAtVariantPosition(vg::ReferenceNode::SharedPtr referenceNodePtr);
		virtual T addReferenceNode(vg::ReferenceNode::SharedPtr referenceNodePtr);

		bool getNextCompoundVariant(Variant::SharedPtr& variant)
		{
			// the first time we call this function we need to get the first variant
			if (!m_next_variant_init)
			{
				m_next_variant_init = true;
				m_variant_list_ptr->getNextVariant(this->m_next_variant);
			}
			variant = this->m_next_variant; // set the variant to be returned
			if (this->m_next_variant == NULL) { return false; } // if the variant is NULL then return false, this indicates we are past the region of interest

			std::string referenceString = variant->getRef()[0];
			std::vector< Variant::SharedPtr > variants;
			Variant::SharedPtr nextVariant;
			bool variantAdded = false;
			position startPosition = variant->getPosition();

			// loop through variants until the variants stop overlapping.
			// As we loop through build a concatenated reference string
			// that represents the entire overlapped variants reference.
			// for those overlapped variants add them to a vector so
			// a "compound variant" can be generated.
			while(m_variant_list_ptr->getNextVariant(nextVariant))
			{
				position variantEndPosition = (variant->getPosition() + referenceString.size() - 1); // subtract 1 because we are counting starting with the position we are on
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
				std::string nextReferenceString = nextVariant->getRef()[0];
				position nextVariantEndPosition = (nextVariant->getPosition() + nextReferenceString.size());
				// if the next variant has reference at a further position then add it to the end of the referenceString
				if (nextVariantEndPosition > variantEndPosition)
				{
					position referenceDelta = (nextVariantEndPosition - variantEndPosition);
					referenceString += nextReferenceString.substr(referenceDelta);
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
			position referenceEndPosition = referenceString.size() + startPosition;
			std::string chrom = variants[0]->getChrom();
			position pos = variants[0]->getPosition();
			std::string id = ".";
			std::string line = chrom + "\t" + std::to_string(pos) + "\t" + id + "\t" + referenceString + "\t";
			// loop over all the variants
			for (auto variantIter = variants.begin(); variantIter != variants.end(); ++variantIter)
			{
				// loop over all the alts in the variants
				for (uint32_t i = 0; i < (*variantIter)->getAlt().size(); ++i)
				{
					std::string altString = (*variantIter)->getAlt()[i];
					// basically we are replacing the variant's reference with the alt within the aggrigated reference (referenceString)
					// it's complicated to explain in words but if you follow the code it isn't too bad
					std::string variantString = referenceString;
					variantString.erase((*variantIter)->getPosition() - startPosition, (*variantIter)->getRef()[0].size());
					variantString.insert((*variantIter)->getPosition() - startPosition, altString);
					line += variantString + ",";
				}
			}
			line.replace(line.size() - 1, 1, "\t\n"); // replace the past comma with a tab
			auto variant = Variant::BuildVariant(line.c_str(), this->m_vcf_parser);
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
