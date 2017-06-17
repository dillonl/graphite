#include "GSSWMapping.h"
#include "core/alignment/BamAlignment.h"

#include <algorithm>
#include <iostream>
#include <thread>
#include <string>

namespace graphite
{
    GSSWMapping::GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr) :
		m_gssw_mapping_ptr(gsswMappingPtr),
        m_offset(m_gssw_mapping_ptr->position),
		m_alignment_ptr(alignmentPtr),
		m_position(0),
		m_mapped(false)
	{
		uint32_t offset = m_gssw_mapping_ptr->position;
		gssw_node_cigar* nc = m_gssw_mapping_ptr->cigar.elements;
		for (int i = 0; i < m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			if (i == 0)
			{
				this->m_position = nc->node->position + m_gssw_mapping_ptr->position;
			}
			auto allelePtr = ((IAllele*)nc->node->data)->getSharedPtr();
			m_allele_ptrs.push_back(allelePtr);
			m_allele_gssw_nodes_map[allelePtr.get()] = nc->node;
		}
	}

	GSSWMapping::~GSSWMapping()
	{
	}

	std::vector< MappingAlignmentInfo::SharedPtr > GSSWMapping::getMappingAlignmentInfoPtrs(IAdjudicator::SharedPtr adjudicatorPtr)
	{
		std::vector< MappingAlignmentInfo::SharedPtr > mappingAlignmentInfoPtrs;
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
        //std::vector< uint32_t > nids;
		for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			int32_t score = 0;
			uint32_t length = 0;
			uint32_t prefixMatch = 0;
			uint32_t suffixMatch = 0;
			auto allelePtr = ((IAllele*)nc->node->data)->getSharedPtr();

            //uint32_t nid = nc->node->id;

			for (int j = 0; j < nc->cigar->length; ++j)
			{
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					score += (adjudicatorPtr->getMatchValue() * nc->cigar->elements[j].length);
					break;
				case 'S':
				case 'X':
					score -= (adjudicatorPtr->getMisMatchValue() * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					score -= adjudicatorPtr->getGapOpenValue();
					score -= (adjudicatorPtr->getGapExtensionValue() * (nc->cigar->elements[j].length -1));
					break;
				default:
					break;
				}
				length += nc->cigar->elements[j].length;
			}
			prefixMatch = length;
			suffixMatch = length;
			score = (score < 0) ? 0 : score; // the floor of the mapping score is 0
			auto mappingAlignmentInfo = std::make_shared< MappingAlignmentInfo >(allelePtr, score, length, prefixMatch, suffixMatch);
			mappingAlignmentInfoPtrs.emplace_back(mappingAlignmentInfo);
		}
		return mappingAlignmentInfoPtrs;
	}

    void GSSWMapping::getGraphPathHeaderAndSequence (std::string& graphPathHeader, std::string& graphPathSequence, position variantPosition)
    {
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;

        std::string nids = "_";
        // Need to alter refOrAlt to account for vcf lines that contain more than one variant.
        int8_t refOrAlt = 0;

        // Loop through nodes to extract node id's and sequences.
        for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
        {
            gssw_node* n = nc->node;
            if (n->id % 2 != 0)
                refOrAlt += 1;
            nids +=  std::to_string(n->id) + ":";
            /*
            if (i != this->m_gssw_mapping_ptr->cigar.length - 1)
                nids += ":";
                */
            graphPathSequence += std::string(n->seq, n->len);
        }

        graphite::BamAlignment::SharedPtr bamAlignmentPtr = std::dynamic_pointer_cast< graphite::BamAlignment >(m_alignment_ptr);
        //graphPathHeader = "RNAME_Header" + std::to_string(bamAlignmentPtr->getRefeqName()) + ":" + std::to_string(variantPosition) + ":" + std::to_string(refOrAlt) + "_" + nids;
        graphPathHeader = bamAlignmentPtr->getAlignmentRegion(bamAlignmentPtr->getRefID()) + ":" + std::to_string(variantPosition) + "->" + "NEW_VARIANT_POSITION" + ":" + std::to_string(refOrAlt) + nids;
    }

    std::string GSSWMapping::getCigarString (IAdjudicator::SharedPtr adjudicatorPtr)
    {
        std::vector< MappingAlignmentInfo::SharedPtr > MappingAlignmentPtrs;
        gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;

        std::string cigar = "";
        for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
        {
            int32_t score = 0;
            uint32_t length = 0;
            uint32_t prefixMatch = 0;
            uint32_t suffixMatch = 0;
            auto allelePtr = ((IAllele*)nc->node->data)->getSharedPtr();

            for (int j = 0; j < nc->cigar->length; ++j)
            {
                cigar += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
            }
        }

        return cigar;
    }

	int GSSWMapping::getMappingScore()
	{
		return this->m_gssw_mapping_ptr->score;
	}

	IAlignment::SharedPtr GSSWMapping::getAlignmentPtr()
	{
		return this->m_alignment_ptr;
	}

	std::vector< IAllele::SharedPtr > GSSWMapping::getAllelePtrs()
	{
		return this->m_allele_ptrs;
	}

	void GSSWMapping::setMapped(bool mapped)
	{
		if (!mapped)
		{
			this->m_allele_incrementor_callback_list.clear();
		}
		this->m_mapped = mapped;
	}

	void GSSWMapping::incrementAlleleCounts()
	{
		if (m_mapped)
		{
			for (auto incrementFunct : this->m_allele_incrementor_callback_list)
			{
				incrementFunct();
			}
		}
	}

	void GSSWMapping::addAlleleCountCallback(std::function< void () > functor)
	{
		this->m_allele_incrementor_callback_list.emplace_back(functor);
	}

	void GSSWMapping::printSimpleMapping()
	{
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
		std::string nodeTypesString = "";
		std::string positions = "";
		std::string alignmentString = std::string(this->m_alignment_ptr->getSequence(), this->m_alignment_ptr->getLength());
		std::string referenceString = "";
		std::string tracebackString = "";
		std::string separator = "||";
		for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			std::string nodeAlignmentSequence = std::string(nc->node->seq, nc->node->len);
			std::string nodeReferenceSequence = std::string(nc->node->ref_seq, nc->node->ref_len);
			positions += std::to_string(nc->node->position) + " ";
			size_t tracebackOffset = 0;
			int refSizeDiff = nc->node->ref_len - nc->node->len;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				uint32_t cigLen = nc->cigar->elements[j].length;
				std::string cigarTracebackString = std::string(nc->node->seq + tracebackOffset, cigLen);
				switch (nc->cigar->elements[j].type)
				{
				case 'S':
				case 'X':
					std::transform(cigarTracebackString.begin(), cigarTracebackString.end(), cigarTracebackString.begin(), ::tolower);
					break;
				case 'M':
				case 'I':
				case 'D':
				default:
					break;
				}
				tracebackString += cigarTracebackString;
			}
			tracebackString += separator;
			nodeTypesString += (nc->node->id % 2 == 0) ? "R" + separator : "V" + separator;
		}
		std::string currentMapping = (this->m_mapped) ? "Mapped" : "Unmapped";
		std::cout << "----------------------------" << std::endl;
		std::cout << "Node Types: " << nodeTypesString.substr(0, nodeTypesString.size() - 2) << std::endl;
		std::cout << "Score:                " << std::to_string(this->m_gssw_mapping_ptr->score) << std::endl;
		std::cout << "Current Mapped State: " << currentMapping << std::endl;
		std::cout << "Positions: " << std::endl;
		std::cout << "Reference: " << referenceString.substr(0, referenceString.size() - 2) << std::endl;
		std::cout << "Traceback: " << tracebackString.substr(0, tracebackString.size() - 2) << std::endl;
		std::cout << "Alignment: " << alignmentString.substr(0, alignmentString.size() - 2) << std::endl;
		std::cout << "----------------------------" << std::endl;
	}

	void GSSWMapping::printMapping()
	{
		position startPosition = 0;
		size_t startSoftClipLength = 0;
		size_t endSoftClipLength = 0;
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
		std::string separator = "||";
		std::string nodeTracebackString = "";
		std::string tracebackString = "";
		std::string cigarString = "";
		std::string alignmentString = this->m_alignment_ptr->getSequence();

		std::vector< position > nodeSeparatorPositions;
		// position referenceStartPosition = 0;
		std::string referenceString = "";
		size_t nodeOffset = 0;
		size_t refOffset = 0;
		size_t nodeSeparatorOffset = 0;

		std::string referenceOffsets;
		std::string tmpAlignmentString;
		std::string tmpTracebackString;
		std::string tmpReferenceString;
		std::string nodeIDString;
		size_t alignmentOffset = 0;
		for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			size_t offset = 0;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				switch (nc->cigar->elements[j].type)
				{
				case 'S':
					if (i == 0 && j == 0)
					{
						startSoftClipLength = nc->cigar->elements[j].length;
						alignmentOffset += startSoftClipLength;
						// tmpAlignmentString += std::string(nc->cigar->elements[j].length, ' ');
					}
					else
					{
						endSoftClipLength = nc->cigar->elements[j].length;
						alignmentString.erase(alignmentString.size() - endSoftClipLength - 1, endSoftClipLength);
					}
					// tmpReferenceString += std::string(nc->node->ref_seq + offset, nc->cigar->elements[j].length);
					break;
				case 'M':
				case 'X':
					tmpAlignmentString += std::string(alignmentString.c_str() + alignmentOffset, nc->cigar->elements[j].length);
					tmpTracebackString += std::string(nc->node->seq + offset, nc->cigar->elements[j].length);
					// tmpReferenceString += std::string(nc->node->ref_seq + offset, nc->cigar->elements[j].length);
					offset += nc->cigar->elements[j].length;
					alignmentOffset += nc->cigar->elements[j].length;
					break;
				case 'I':
					tmpAlignmentString += std::string(alignmentString.c_str() + alignmentOffset, nc->cigar->elements[j].length);
					tmpTracebackString += std::string(nc->cigar->elements[j].length, '-');
					// tmpReferenceString += std::string(nc->node->ref_seq + offset, nc->cigar->elements[j].length);
					offset += nc->cigar->elements[j].length;
					alignmentOffset += nc->cigar->elements[j].length;
					break;
				case 'D':
					tmpAlignmentString += std::string(alignmentString.c_str() + alignmentOffset, nc->cigar->elements[j].length);
					// tmpReferenceString += std::string(nc->cigar->elements[j].length, '-');
					tmpTracebackString += std::string(nc->node->seq + offset, nc->cigar->elements[j].length);
					offset += nc->cigar->elements[j].length;
					alignmentOffset += nc->cigar->elements[j].length;
					break;
				}
				cigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
			}
			if (startPosition == 0)
			{
				startPosition = nc->node->position;
			}
			std::string nodeTypeString = (nc->node->id % 2 == 0) ? "R" : "V";
			nodeIDString += std::to_string(nc->node->id) + "||";

			cigarString +=  separator;
			tmpAlignmentString += separator;
			tmpTracebackString += separator;
			tmpReferenceString += std::string(nc->node->ref_seq, nc->node->ref_len) + separator;

			referenceString += nc->node->ref_seq;
			tracebackString += nc->node->seq;
			nodeSeparatorPositions.emplace_back(nc->node->position + nodeSeparatorOffset);

			if (nc->node->ref_len > nc->node->len) { nodeSeparatorOffset = nc->node->ref_len - nc->node->len; tracebackString += std::string(nc->node->ref_len - nc->node->len, ' '); }
			else if (nc->node->ref_len < nc->node->len) { nodeSeparatorOffset = nc->node->len - nc->node->ref_len; referenceString += std::string(nc->node->len - nc->node->ref_len, ' '); }
			if (i < this->m_gssw_mapping_ptr->cigar.length && i != this->m_gssw_mapping_ptr->cigar.length - 1)
			{
				referenceString += separator;
				tracebackString += separator;
			}

			nodeTracebackString += nodeTypeString + separator;
			nodeOffset += nc->node->len - 1;
			refOffset += nc->node->ref_len;
		}
		nodeTracebackString = (nodeTracebackString.size() > 2) ? nodeTracebackString.substr(0, nodeTracebackString.size() - 2) : nodeTracebackString;
		cigarString = (cigarString.size() > 2) ? cigarString.substr(0, cigarString.size() - 2) : cigarString;
		tmpAlignmentString = (tmpAlignmentString.size() > 2) ? tmpAlignmentString.substr(0, tmpAlignmentString.size() - 2) : tmpAlignmentString;
		tmpTracebackString = (tmpTracebackString.size() > 2) ? tmpTracebackString.substr(0, tmpTracebackString.size() - 2) : tmpTracebackString;
		tmpReferenceString = (tmpReferenceString.size() > 2) ? tmpReferenceString.substr(0, tmpReferenceString.size() - 2) : tmpReferenceString;

		alignmentString = std::string(this->m_gssw_mapping_ptr->position, ' ') + alignmentString.substr(startSoftClipLength);

		size_t sepPos = 0;
		while ((sepPos = tracebackString.find(separator, sepPos + 1)) != std::string::npos && sepPos < alignmentString.size())
		{
			alignmentString.insert(sepPos, separator);
		}

		std::string isMapped = (this->m_alignment_ptr->isMapped()) ? "Mapped" : "Unmapped";
		std::string mate = (this->m_alignment_ptr->isFirstMate()) ? "First" : "Second";
		std::string reverseStrand = (this->m_alignment_ptr->isReverseStrand()) ? "Reverse" : "Normal";
		std::string currentMapping = (this->m_mapped) ? "Mapped" : "Unmapped";

		std::string reportString = "";
		std::string eol = "\r\n";
		reportString += "---------------------------------------------------------" + eol;
		reportString += "Score:                " + std::to_string(this->m_gssw_mapping_ptr->score) + eol;
		reportString += "Node Cigar String:    " + cigarString + eol;
		reportString += "Node Type Traceback:  " + nodeTracebackString + eol;
		reportString += "Node ID Traceback:    " + nodeIDString + eol;
		reportString += "Reference:            " + referenceString + eol;
		// reportString += "Reference (new):      " + tmpReferenceString + eol;
		// reportString += "Alignment (new):      " + tmpAlignmentString + eol;
		// reportString += "Traceback (new):      " + tmpTracebackString + eol;
		reportString += "Node Traceback:       " + tracebackString + eol;
		reportString += "Alignment:            " + alignmentString + eol;
		reportString += "Alignment(RAW):       " + std::string(this->m_alignment_ptr->getSequence(), this->m_alignment_ptr->getLength()) + eol;
		reportString += "Alignment Name:       " + this->m_alignment_ptr->getID() + eol;
		reportString += "Alignment Position:   " + std::to_string(this->m_alignment_ptr->getPosition()) + eol;
		reportString += "Graph Position:       " + std::to_string(startPosition) + eol;
		reportString += "Mapping ID:           " + std::to_string(this->m_id) + eol;
		reportString += "Current Mapped State: " + currentMapping + eol;
		reportString += "Prev. Mapped State:   " + isMapped + eol;
		reportString += "Mate Order:           " + mate + eol;
		reportString += "Strand:               " + reverseStrand + eol;
		reportString += "Original Map Quality: " + std::to_string(this->m_alignment_ptr->getOriginalMapQuality()) + eol;
		reportString += "---------------------------------------------------------" + eol;

		std::cout << reportString;
	}

}
