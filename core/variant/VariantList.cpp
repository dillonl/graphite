#include "VariantList.h"
#include "Variant.h"

namespace gwiz
{
	VariantList::VariantList(const std::vector< IVariant::SharedPtr >& variantPtrs) :
		m_variant_ptrs(variantPtrs),
		m_current_index(0),
		m_next_variant_init(false)
	{
	}

	VariantList::~VariantList()
	{
	}

	bool VariantList::getNextVariant(IVariant::SharedPtr& variantPtr)
	{
		if (this->m_current_index < this->m_variant_ptrs.size())
		{
			variantPtr = this->m_variant_ptrs[this->m_current_index++];
			return true;
		}
		else
		{
			variantPtr = nullptr;
			return false;
		}
	}

	size_t VariantList::getCount()
	{
		return this->m_variant_ptrs.size();
	}

	void VariantList::sort()
	{
		std::sort(this->m_variant_ptrs.begin(), this->m_variant_ptrs.end(),
				  [] (const IVariant::SharedPtr a, const IVariant::SharedPtr b)
				  {
					  return a->getPosition() < b->getPosition();
				  }
			);
	}

	void VariantList::normalizeOverlappingVariants()
	{
		std::vector< IVariant::SharedPtr > variantPtrs;
		IVariant::SharedPtr variantPtr;
		while (getNextCompoundVariant(variantPtr))
		{
			variantPtrs.emplace_back(variantPtr);
		}
		this->m_variant_ptrs = variantPtrs;
	}

	bool VariantList::hasVariants()
	{
		return getCount() > 0;
	}

	void VariantList::addVariants(IVariantList::SharedPtr variantListPtr)
	{
		throw "VariantList::addVariants not implemented";
	}

	void VariantList::printHeader(std::ostream& out)
	{
		out << "##fileformat=VCFv4.2" << std::endl;
		out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">" << std::endl;
		out << "##INFO=<ID=TC,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">" << std::endl;
		out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth including non-counted low quality reads\">" << std::endl;
		out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
		out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" << std::endl;
	}

	void VariantList::printToVCF(std::ostream& out)
	{
		printHeader(out);
		for(const auto variantPtr : this->m_variant_ptrs)
		{
			variantPtr->printVariant(out);
		}
	}

	bool VariantList::getNextCompoundVariant(IVariant::SharedPtr& variant)
	{
		// the first time we call this function we need to get the first variant
		if (!m_next_variant_init)
		{
			m_next_variant_init = true;
			getNextVariant(this->m_next_variant);
		}
		variant = this->m_next_variant; // set the variant to be returned
		if (this->m_next_variant == NULL) { return false; } // if the variant is NULL then return false, this indicates we are past the region of interest

		std::string referenceString = variant->getRefAllelePtr()->getSequenceString();
		std::vector< IVariant::SharedPtr > variants;
		IVariant::SharedPtr nextVariant;
		bool variantAdded = false;
		position startPosition = variant->getPosition();
		std::string variantChrom = variant->getChrom();
		position compoundStartPosition = startPosition;
		position variantEndPosition = (startPosition + referenceString.size() - 1); // subtract 1 because we are counting starting with the position we are on

		// loop through variants until the variants stop overlapping.
		// As we loop through build a concatenated reference string
		// that represents the entire overlapped variants reference.
		// for those overlapped variants add them to a vector so
		// a "compound variant" can be generated.
		while (getNextVariant(nextVariant))
		{
			if (variantEndPosition < nextVariant->getPosition() || strcmp(variantChrom.c_str(), nextVariant->getChrom().c_str()) != 0)
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
			std::string nextReferenceString = nextVariant->getRefAllelePtr()->getSequence();
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
			variantChrom = nextVariant->getChrom();
		}

		if (!variants.empty())
		{
			variant = buildCompoundVariant(startPosition, referenceString, variants);
		}
		this->m_next_variant = nextVariant; // set the next variant
		return true;
	}

	IVariant::SharedPtr VariantList::buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< IVariant::SharedPtr >& variants)
	{
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
			for (auto altAllelePtr : variantPtr->getAltAllelePtrs())
			{
				std::string altString = std::string(altAllelePtr->getSequence());
				// basically we are replacing the variant's reference with the alt within the aggrigated reference (referenceString)
				// it's complicated to explain in words but if you follow the code it isn't too bad
				std::string variantString = referenceString;
				variantString.erase(variantPtr->getPosition() - startPosition, variantPtr->getRefAllelePtr()->getSequenceString().size());
				variantString.insert(variantPtr->getPosition() - startPosition, altString);
				line += variantString + ",";
			}
		}
		line.replace(line.size() - 1, 1, "\t"); // replace the past comma with a tab
		std::string quality = "0";
		std::string filter = "PASS";
		std::string info = "C=TRUE";
		line += quality + "\t" + filter + "\t" + info + "\t.\n";
		auto variant = Variant::BuildVariant(line.c_str(), this->m_vcf_parser);
		return variant;
	}

	VariantList::SharedPtr VariantList::getVariantsInRegion(Region::SharedPtr regionPtr)
	{
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();


		auto lowerBound = std::lower_bound(this->m_variant_ptrs.begin(), this->m_variant_ptrs.end(), nullptr, [startPosition](const IVariant::SharedPtr& variantPtr, const IVariant::SharedPtr& ignore) {
				return startPosition > variantPtr->getPosition();
			});
		auto upperBound = std::upper_bound(this->m_variant_ptrs.begin(), this->m_variant_ptrs.end(), nullptr, [endPosition](const IVariant::SharedPtr& ignore, const IVariant::SharedPtr& variantPtr) {
				return variantPtr->getPosition() > endPosition;
			});

		auto variantPtrs = std::vector< IVariant::SharedPtr >(lowerBound, upperBound);
		return std::make_shared< VariantList >(variantPtrs);
	}
}
