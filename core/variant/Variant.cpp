#include "Variant.h"
#include "core/alignment/IAlignment.h"
#include "core/util/Utility.h"

#include <unordered_map>
#include <sstream>
// #include <regex>

namespace graphite
{
	Variant::Variant(position pos, const std::string& chrom, const std::string& id, const std::string& quality, const std::string& filter, IAllele::SharedPtr refAllelePtr, std::vector< IAllele::SharedPtr > altAllelePtrs, uint32_t readLength) : m_position(pos), m_chrom(chrom), m_id(id), m_qual(quality), m_filter(filter), m_skip(false), m_read_length(readLength), m_reference_size(0), m_overlap(10)
	{
		this->m_ref_allele_ptr = refAllelePtr;
		this->m_alt_allele_ptrs = altAllelePtrs;

		this->m_all_allele_ptrs.reserve(this->m_alt_allele_ptrs.size() + 1);
		this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
		this->m_all_allele_ptrs.insert(this->m_all_allele_ptrs.end(), this->m_alt_allele_ptrs.begin(), this->m_alt_allele_ptrs.end());
	}

	Variant::Variant() :
		m_skip(false), m_overlap(10)
	{
	}

	Variant::~Variant()
	{
	}

	Variant::SharedPtr Variant::BuildVariant(const std::string& vcfLine, IReference::SharedPtr referencePtr, uint32_t readLength)
	{
		std::string fields;
		auto variantPtr = std::make_shared< Variant >();
		std::string ref;
		std::string alts;
		std::vector< std::string > vcfComponents;

		split(vcfLine, '\t', vcfComponents);

		if (vcfComponents.size() < 7)
		{
			std::cout << "vcf line is incorrectly formated" << std::endl;
		}

		variantPtr->m_chrom = vcfComponents[0];
		variantPtr->m_position = stoul(vcfComponents[1]);
		variantPtr->m_id = vcfComponents[2];
		ref = vcfComponents[3];
		alts = vcfComponents[4];

		{
			static std::mutex l;
			std::lock_guard< std::mutex > lg(l);
			std::cout << "working on: " << variantPtr->m_chrom << " " << variantPtr->m_position << std::endl;
		}

		variantPtr->m_qual = vcfComponents[5];
		variantPtr->m_filter = vcfComponents[6];
		fields = vcfComponents[7];

		variantPtr->m_info_fields.clear();
		variantPtr->setUnorderedMapKeyValue(fields);
		variantPtr->setAlleles(referencePtr, ref, alts, readLength);

		variantPtr->m_line = std::string(vcfLine.c_str());

		return variantPtr;
	}

	void Variant::setUnorderedMapKeyValue(const std::string& rawString)
	{
		m_info_fields.clear();
		std::vector< std::string > infoFieldKeyValueSplit;
		split(rawString, ';', infoFieldKeyValueSplit);
		for (const auto& infoField : infoFieldKeyValueSplit)
		{
			std::vector< std::string > infoFieldSplit;
			split(infoField, '=', infoFieldSplit);
			if (infoFieldSplit.size() == 2)
			{
				m_info_fields[infoFieldSplit[0]] = infoFieldSplit[1];
			}
		}
	}

	int getSVLengthFromInfo(std::unordered_map< std::string, std::string >& infoFields, position pos)
	{
		int svLength = 0;
		if (infoFields.find("SVLEN") != infoFields.end())
		{
			svLength = stoi(infoFields["SVLEN"]);
		}
		else if (infoFields.find("END") != infoFields.end())
		{
			svLength = stoi(infoFields["END"]) - pos;
		}
		else if (infoFields.find("SEQ") != infoFields.end())
		{
			svLength = infoFields["SEQ"].size();
		}
		if (svLength < 0) { svLength = svLength * -1; }
		return svLength;
	}

	std::string getInsertionSequenceFromInfo(std::unordered_map< std::string, std::string >& infoFields)
	{
		auto seqIter = infoFields.find("SEQ");
		if (seqIter == infoFields.end())
		{
			return "";
		}
		else
		{
			return seqIter->second;
		}
	}

	bool isStandardAlt(const std::string& alt)
	{
		std::unordered_set< char > validSequence = {'A', 'C', 'G', 'T', ','};
		for (auto c : alt)
		{
			if (validSequence.find(c) == validSequence.end())
			{
				return false;
			}
		}
		return true;
	}

	bool isSymbolicAlt(const std::string& alt)
	{
		std::unordered_set< char > symbolicSequence = {'<', '>'};
		for (auto c : alt)
		{
			if (symbolicSequence.find(c) != symbolicSequence.end())
			{
				return true;
			}
		}
		return false;
	}

	std::string getTruncatedSequence(const char* sequence, uint32_t svLength, uint32_t readLength)
	{
		std::cout << std::string(sequence, readLength) + std::string(readLength, 'N') + std::string(sequence + (svLength - readLength), readLength) << std::endl;
		return std::string(sequence, readLength) + std::string(readLength, 'N') + std::string(sequence + (svLength - readLength), readLength);
	}

	void Variant::setAsDeletion(Reference::SharedPtr referencePtr, int svLength, uint32_t readLength)
	{
		auto maxSize = readLength * 3;
		// auto maxSize = std::numeric_limits< uint32_t >::max();
		std::string ref;
		if (svLength > maxSize) // if the variant is too large then create a truncated reference allele with Ns seperating the two breakpoints
		{
			ref = referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, m_position, m_position + readLength, Region::BASED::ONE)) + std::string(readLength, 'N') + referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, (m_position + 1) + svLength - readLength, m_position + svLength, Region::BASED::ONE));
			m_reference_size = svLength;
		}
		else
		{
			ref = referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, m_position, (m_position + svLength), Region::BASED::ONE));
			m_reference_size = ref.size();
		}
		std::vector< std::string > alts;
		alts.emplace_back(std::string(1, ref[0]));
		setRefAndAltAlleles(ref, alts);
	}

	void Variant::setAsDuplication(Reference::SharedPtr referencePtr, int svLength, uint32_t readLength)
	{
		auto maxSize = readLength * 3;
		// auto maxSize = std::numeric_limits< uint32_t >::max();
		const char* reference = referencePtr->getSequence() + (m_position - referencePtr->getRegion()->getStartPosition());
		std::string ref;
		std::string alt;
		std::vector< std::string > alts;
		if (svLength > maxSize)
		{
			std::string sequenceA = std::string(reference + 1, readLength);
			std::string intermediateSequence(readLength, 'N');
			std::string sequenceB = std::string((reference + 1) + (svLength - readLength), readLength);
			ref = std::string(1, reference[0]) + sequenceA + intermediateSequence + sequenceB;
			alt = std::string(1, reference[0]) + sequenceA + intermediateSequence + sequenceB + sequenceA + intermediateSequence + sequenceB;

			m_reference_size = svLength;
		}
		else
		{
			ref = std::string(reference, svLength + 1);
			alt = ref + std::string(reference + 1, svLength); // the duplicated region does not contain the first bp of the reported ref sequence

			m_reference_size = ref.size();
		}
		alts.emplace_back(alt);
		setRefAndAltAlleles(ref, alts);
	}

	void Variant::setAsInversion(Reference::SharedPtr referencePtr, int svLength, uint32_t readLength)
	{
		auto maxSize = readLength * 3;
		// auto maxSize = std::numeric_limits< uint32_t >::max();
		const char* reference = referencePtr->getSequence() + (m_position - referencePtr->getRegion()->getStartPosition());
		std::string ref;
		std::string alt;
		if (svLength > maxSize)
		{
			std::string intermediateSequence(readLength, 'N');
			std::string anchorBase = referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, m_position, (m_position), Region::BASED::ONE));
			std::string sequenceA = referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, (m_position + 1), (m_position + readLength), Region::BASED::ONE));
			std::string sequenceB = referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, ((m_position + 1) + svLength - readLength), (m_position + svLength), Region::BASED::ONE));
			ref = anchorBase + sequenceA + intermediateSequence + sequenceB;
			std::reverse(sequenceA.begin(), sequenceA.end());
			std::reverse(sequenceB.begin(), sequenceB.end());
			alt = anchorBase + sequenceB + intermediateSequence + sequenceA;
			m_reference_size = svLength;
		}
		else
		{
			ref = referencePtr->getSequenceFromRegion(std::make_shared< Region >(this->m_chrom, m_position, (m_position + svLength), Region::BASED::ONE));
			alt = std::string(ref.c_str() + 1, svLength);
			std::reverse(alt.begin(), alt.end());
			alt = std::string(1, ref.c_str()[0]) + alt; // add on the anchor base

			ref.size();
		}
		std::vector< std::string > alts;
		alts.emplace_back(alt);
		setRefAndAltAlleles(ref, alts);
	}

	void Variant::setAsInsertion(const std::string& ref, const std::string& alt, uint32_t readLength)
	{
		auto maxSize = readLength * 3;
		// auto maxSize = std::numeric_limits< uint32_t >::max();
		std::vector< std::string > alts;
		if (alt.size() > maxSize)
		{
			// std::string altSequence = ref + getTruncatedSequence(alt.c_str(), alt.size(), readLength);
			std::string altSequence = getTruncatedSequence(alt.c_str(), alt.size(), readLength);
			alts.emplace_back(altSequence);
		}
		else
		{
			alts.emplace_back(alt);
		}
		setRefAndAltAlleles(ref, alts);
		m_reference_size = ref.size();
	}

	void Variant::setAsStandardAlt(const std::string& ref, const std::string& alt, uint32_t readLength)
	{
		m_reference_size = ref.size();
		std::string tmpRef;
		m_variant_size = ref.size();
		if (ref.size() > (readLength * 3))
		{
			m_is_sv = true;
			tmpRef = ref.c_str()[0] + getTruncatedSequence(ref.c_str() + 1, ref.size() - 1, readLength);
		}
		else
		{
			tmpRef = ref;
		}

		std::vector< std::string > tmpAlts;
		std::vector< std::string > alts;
		if (alt.find(",") != std::string::npos)
		{
			split(alt, ',', tmpAlts);
		}
		else
		{
			tmpAlts.emplace_back(alt);
		}
		// for each alt allele make sure it's not larger than the maxSize, if it is then truncate it and add the truncated allele
		for (auto& tmpAlt : tmpAlts)
		{
			if (tmpAlt.size() > m_variant_size)
			{
				m_variant_size = tmpAlt.size();
			}
			if (tmpAlt.size() > (readLength * 3))
			{
				m_is_sv = true;
				alts.emplace_back(getTruncatedSequence(tmpAlt.c_str(), tmpAlt.size(), readLength));
			}
			else
			{
				alts.emplace_back(tmpAlt);
			}
		}

		setRefAndAltAlleles(tmpRef, alts);
	}

	void Variant::setAlleles(Reference::SharedPtr referencePtr, const std::string& vcfReferenceString, const std::string& alt, uint32_t readLength)
	{
		m_is_sv = false;
		bool shouldSkip = true;
		int svLength = getSVLengthFromInfo(m_info_fields, m_position);
		if (isStandardAlt(alt))
		{
			setAsStandardAlt(vcfReferenceString, alt, readLength);
			shouldSkip = false;
			svLength = vcfReferenceString.size();
			m_variant_size = svLength;
		}
		else if (isSymbolicAlt(alt) && svLength > 0)
		{
			m_is_sv = (svLength > (readLength * 3));
			m_variant_size = svLength;
			if (alt.compare("<DEL>") == 0)
			{
				setAsDeletion(referencePtr, svLength, readLength);
				shouldSkip = false;
			}
			else if (alt.compare("<DUP>") == 0)
			{
				setAsDuplication(referencePtr, svLength, readLength);
				shouldSkip = false;
			}
			else if (alt.compare("<INV>") == 0)
			{
				setAsInversion(referencePtr, svLength, readLength);
				shouldSkip = false;
			}
			else if (alt.compare("<INS>") == 0)
			{
				std::string insertionSequence  = getInsertionSequenceFromInfo(m_info_fields);
				if (insertionSequence.size() > 0)
				{
					setAsInsertion(vcfReferenceString, insertionSequence, readLength);
					shouldSkip = false;
				}

			}
		}
		if (shouldSkip)
		{
			setSkip(true);
		}
		else
		{
			if (m_is_sv)
			{
				position startPosition1 = (this->m_position - (readLength + m_overlap) <= 0) ? 0 : (this->m_position - (readLength + m_overlap));
				position endPosition1 = this->m_position + (readLength + m_overlap);
				position startPosition2 = (((this->m_position + 1) + svLength - (readLength + m_overlap)) <= 0) ? 0 : ((this->m_position + 1) + svLength - (readLength + m_overlap));
				position endPosition2 = (this->m_position + 1 + svLength + readLength + m_overlap);
				m_region_ptrs.emplace_back(std::make_shared< Region >(m_chrom, startPosition1, endPosition1, Region::BASED::ONE));
				m_region_ptrs.emplace_back(std::make_shared< Region >(m_chrom, startPosition2, endPosition2, Region::BASED::ONE));
			}
			else
			{
				position startPosition1 = (this->m_position - (readLength + m_overlap) <= 0) ? 0 : (this->m_position - (readLength + m_overlap));
				position endPosition1 = (this->m_position + svLength + 1 + readLength + m_overlap);
				m_region_ptrs.emplace_back(std::make_shared< Region >(m_chrom, startPosition1, endPosition1, Region::BASED::ONE));
			}
		}
	}

	uint32_t Variant::getAllelePrefixOverlapMaxCount(IAllele::SharedPtr allelePtr)
	{
		auto iter = this->m_allele_prefix_max_overlap_map.find(allelePtr);
		return (iter != this->m_allele_prefix_max_overlap_map.end()) ? iter->second : 0;
	}

	uint32_t Variant::getAlleleSuffixOverlapMaxCount(IAllele::SharedPtr allelePtr)
	{
		auto iter = this->m_allele_suffix_max_overlap_map.find(allelePtr);
		return (iter != this->m_allele_suffix_max_overlap_map.end()) ? iter->second : 0;
	}

	void Variant::setAlleleOverlapMaxCountIfGreaterThan(IAllele::SharedPtr allelePtr, std::unordered_map< IAllele::SharedPtr, uint32_t >& alleleOverlapCountMap, uint32_t overlapCount)
	{
		auto iter = alleleOverlapCountMap.find(allelePtr);
		if (iter == alleleOverlapCountMap.end() || iter->second < overlapCount)
		{
			alleleOverlapCountMap[allelePtr] = overlapCount;
		}
	}

	/*
	 * This is where allele's variant weakptr is set.
	 * Also, it is with respect to this variant that
	 * the allele prefix and suffix match is calculated.
	 */
	void Variant::processOverlappingAlleles()
	{
		std::weak_ptr< IVariant > weakPtr = shared_from_this();
		uint32_t alleleCount = 1;
		m_max_prefix_match_length = 0;
		m_max_suffix_match_length = 0;
		for (auto& allelePtr1 : this->m_all_allele_ptrs)
		{
			allelePtr1->setVariantWPtr(weakPtr);
			for (uint32_t i = alleleCount; i < this->m_all_allele_ptrs.size(); ++i)
			{
				auto allelePtr2 = this->m_all_allele_ptrs[i];
				auto tmpMaxPrefixCount = allelePtr1->getCommonPrefixSize(allelePtr2);
				auto tmpMaxSuffixCount = allelePtr1->getCommonSuffixSize(allelePtr2);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr1, this->m_allele_prefix_max_overlap_map, tmpMaxPrefixCount);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr2, this->m_allele_prefix_max_overlap_map, tmpMaxPrefixCount);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr1, this->m_allele_suffix_max_overlap_map, tmpMaxSuffixCount);
				setAlleleOverlapMaxCountIfGreaterThan(allelePtr2, this->m_allele_suffix_max_overlap_map, tmpMaxSuffixCount);
			}
			++alleleCount;
		}
	}

	void Variant::setRefAndAltAlleles(const std::string& ref, const std::vector< std::string >& alts)
	{
		int id = 0;
		this->m_all_allele_ptrs.clear();
		std::string upperRef = ref;
		std::transform(upperRef.begin(), upperRef.begin(), upperRef.end(), ::toupper);
		this->m_ref_allele_ptr = std::make_shared< Allele >(upperRef);
		this->m_ref_allele_ptr->setID(id);
		this->m_all_allele_ptrs.reserve(alts.size() + 1);
		this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
		this->m_alt_allele_ptrs.clear();
		id += 1;
		for (const auto& alt : alts)
		{
			std::string upperAlt = alt;
			std::transform(upperAlt.begin(), upperAlt.begin(), upperAlt.end(), ::toupper);
			auto altAllelePtr = std::make_shared< Allele >(upperAlt);
			altAllelePtr->setID(id);
			this->m_alt_allele_ptrs.emplace_back(altAllelePtr);
			this->m_all_allele_ptrs.emplace_back(altAllelePtr);
			id += 2;
		}
	}

	bool Variant::shouldSkip() { return this->m_skip; }

	void Variant::setSkip(bool skip)
	{
		m_skip = skip;
		std::string ref = "";
		std::vector< std::string > alts;
		setRefAndAltAlleles(ref, alts);
	}

	bool Variant::doesOverlap(IVariant::SharedPtr variantPtr)
	{
		auto variantRegionPtrs =  variantPtr->getRegions();
		for (auto thisRegionPtr : m_region_ptrs)
		{
			auto thisStartPosition = thisRegionPtr->getStartPosition();
			auto thisEndPosition = thisRegionPtr->getEndPosition();
			for (auto variantRegionPtr : variantRegionPtrs)
			{
				auto variantStartPosition = variantRegionPtr->getStartPosition();
				auto variantEndPosition = variantRegionPtr->getEndPosition();
				if ((thisStartPosition <= variantStartPosition && variantStartPosition <= thisEndPosition) ||
					(variantStartPosition <= thisStartPosition && thisStartPosition <= variantEndPosition))
				{
					return true;
				}
			}
		}
		return false;
	}

	uint32_t Variant::getReferenceSize()
	{
		return m_reference_size;
	}

	void Variant::addRegion(Region::SharedPtr regionPtr)
	{
		m_region_ptrs.emplace_back(regionPtr);
	}

	std::string Variant::getInfoFieldsString()
	{
		std::string infoFields = "";
		for (auto infoField : this->m_info_fields)
		{
			std::string prefix = (infoFields.size() > 0) ? ";" : "";
			infoFields += prefix + infoField.first + "=" + infoField.second;
		}
		return (this->m_info_fields.size() > 0) ? infoFields : ".";
	}

	std::string Variant::getVariantLine(IHeader::SharedPtr headerPtr)
	{
		std::vector< std::string > lineSplit;
		split(this->m_line, '\t', lineSplit);
		std::string line = "";

		std::string samplePadding = "";
		bool samplePaddingSet = false;
		for (auto i = lineSplit.size(); i < 9 + headerPtr->getColumnNames().size(); ++i) // add blank sample columns
		{
			lineSplit.emplace_back("");

			if (!samplePaddingSet)
			{
				samplePaddingSet = true;
				auto formatIdx = headerPtr->getColumnPosition("FORMAT");
				auto formatField = lineSplit[formatIdx];
				auto numFields = std::count(formatField.begin(), formatField.end(), ':');
				if (formatField.size() > 0)
				{
					samplePadding += ".";
					for (auto n = 0; n < numFields; ++n) { samplePadding += ":."; }
				}
			}
		}

	    auto formatColumnIdx = headerPtr->getColumnPosition("FORMAT");
		uint32_t i = 0;
		for (i = 0; i < 9; ++i)
		{
			line += (i > 0) ? "\t" : "";
			if (i == formatColumnIdx) // append the format info to the format columns
			{
				if (lineSplit.size() >= i && !lineSplit[i].empty())
				{
					line += lineSplit[i] + ":";
				}
				line += getFormatString();
			}
			else if (lineSplit.size() >= i)
			{
				line += lineSplit[i];
			}
		}

		for (; i < headerPtr->getColumnNames().size(); ++i)
		{
			auto columnName = headerPtr->getColumnNames()[i];
			line += "\t";
			if (lineSplit[i].empty())
			{
				line += (samplePadding.size() > 0) ? samplePadding + ":" : "";
			}
			else
			{
				line += lineSplit[i] + ":";
			}
			auto sampleCounts = (headerPtr->isActiveSampleColumnName(columnName)) ?  getSampleCounts(columnName) : getBlankSampleCounts();
			line += sampleCounts;
		}

		return line + "\n";
	}

	std::string Variant::getFormatString()
	{
		std::string formatString = "";
		for (auto i = 0; i < AllAlleleCountTypes.size(); ++i)
		{
			auto alleleCountType = AllAlleleCountTypes[i];
			std::string suffix = (i < AllAlleleCountTypes.size() - 1) ? ":" : "";
			std::string alleleTypeCountString = AlleleCountTypeToShortString(alleleCountType);
			formatString += "DP_" + alleleTypeCountString + ":DP4_" + alleleTypeCountString + suffix;
		}
		return formatString;
	}

	std::string Variant::getBlankSampleCounts()
	{
		std::string alleleCountString = "";
	    for (auto i = 0; i < AllAlleleCountTypes.size(); ++i)
		{
			alleleCountString += ".:.:";
		}
		alleleCountString += ".";
		return alleleCountString;
	}

	std::string Variant::getSampleCounts(const std::string& sampleName)
	{
		std::string alleleCountString = "";
		for (auto i = 0; i < AllAlleleCountTypes.size(); ++i)
		{
			auto alleleCountType = AllAlleleCountTypes[i];
			std::string suffix = (i < AllAlleleCountTypes.size() - 1) ? ":" : "";
			std::string alleleTypeCountString = AlleleCountTypeToString(alleleCountType);
			uint32_t totalCount = 0;
			std::string sampleString = "";
			for (size_t j = 0; j < m_all_allele_ptrs.size(); ++j)
			{
				auto allelePtr = m_all_allele_ptrs[j];
				uint32_t forwardCount = allelePtr->getForwardCount(sampleName, alleleCountType);
				uint32_t reverseCount = allelePtr->getReverseCount(sampleName, alleleCountType);
				std::string prefix = (j == 0) ? "" : ",";

				if (m_skip)
				{
					sampleString += prefix + ".,.";
				}
				else
				{
					sampleString += prefix + std::to_string(forwardCount) + "," + std::to_string(reverseCount);
				}

				totalCount += (forwardCount + reverseCount);
			}
			std::string totalCountString = (m_skip) ? "." : std::to_string(totalCount);
			// alleleCountString += "DP<" + alleleTypeCountString + ">=" + totalCountString + ";DP4<" + alleleTypeCountString + ">=" + sampleString + ";";
			alleleCountString += totalCountString + ":" + sampleString + suffix;
		}
		return alleleCountString;
	}

}
