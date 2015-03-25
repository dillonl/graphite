#include "AlignmentReport.h"
#include "core/genotyper/GenotyperAllele.hpp"

namespace gwiz
{
namespace gssw
{
	AlignmentReport::AlignmentReport(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignment::SharedPtr alignmentPtr, std::shared_ptr< gssw_graph_mapping > graphMappingPtr, position graphStartPosition) :
		m_reference_ptr(referencePtr),
		m_variant_list_ptr(variantListPtr),
		m_alignment_ptr(alignmentPtr),
		m_graph_mapping_ptr(graphMappingPtr),
		m_graph_start_position(graphStartPosition)
	{
	}

	AlignmentReport::~AlignmentReport()
	{
	}

	std::string AlignmentReport::toString()
	{
		position startPosition = 0;
		size_t startSoftClipLength = 0;
		size_t endSoftClipLength = 0;
		gssw_node_cigar* nc = this->m_graph_mapping_ptr->cigar.elements;
		std::string separator = "||";
		std::string nodeTracebackString = "";
		std::string tracebackString = "";
		std::string cigarString = "";
		std::string alignmentString = this->m_alignment_ptr->getSequence();
		std::vector< position > nodeSeparatorPositions;
		size_t nodeOffset = 0;
		for (int i = 0; i < this->m_graph_mapping_ptr->cigar.length; ++i, ++nc)
		{
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				if (nc->cigar->elements[j].type == 'S') // softclipping, set softclip variables
				{
					if (i == 0 && j == 0)
					{
						startSoftClipLength = nc->cigar->elements[j].length;
					}
					else
					{
						endSoftClipLength = nc->cigar->elements[j].length;
					}
				}
				cigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
			}
			cigarString +=  separator;
			nodeSeparatorPositions.emplace_back(nc->node->position);
			// nodeSeparatorIndices.emplace_back(i + nodeOffset);
			if (startPosition == 0) { startPosition = nc->node->position; }
			tracebackString += nc->node->seq;
			nodeTracebackString += std::string(GenotyperAllele::TypeToString(static_cast< GenotyperAllele::Type >((long)nc->node->data))) + separator;
			nodeOffset += nc->node->len - 1;
		}
		nodeTracebackString = (nodeTracebackString.size() > 2) ? nodeTracebackString.substr(0, nodeTracebackString.size() - 2) : nodeTracebackString;
		cigarString = (cigarString.size() > 2) ? cigarString.substr(0, cigarString.size() - 2) : cigarString;

		position referenceStartPosition = startPosition - this->m_reference_ptr->getRegion()->getStartPosition();
		std::string referenceString = std::string((this->m_reference_ptr->getSequence() + referenceStartPosition), tracebackString.size());
		alignmentString = std::string(this->m_graph_mapping_ptr->position, ' ') + alignmentString.substr(startSoftClipLength);

		for (int i = nodeSeparatorPositions.size() - 1; i > 0; --i)
		{
			// size_t index = nodeSeparatorIndices.at(i);
			size_t index = nodeSeparatorPositions.at(i) - startPosition;
			if (tracebackString.size() > index) { tracebackString.insert(index, separator); }
			if (referenceString.size() > index) { referenceString.insert(index, separator); }
			if (alignmentString.size() > index) { alignmentString.insert(index, separator); }
		}

		// if (true) { return ""; }

		std::string reportString = "";
		std::string eol = "\r\n";
		reportString += "---------------------------------------------------------" + eol;
		reportString += "Score:                " + std::to_string(this->m_graph_mapping_ptr->score) + eol;
		reportString += "Node Cigar String:    " + cigarString + eol;
		reportString += "Node Type Traceback:  " + nodeTracebackString + eol;
		reportString += "Reference:            " + referenceString + eol;
		reportString += "Node Traceback:       " + tracebackString + eol;
		reportString += "Alignment:            " + alignmentString + eol;
		reportString += "Alignment Position:   " + std::to_string(this->m_alignment_ptr->getPosition()) + eol;
		reportString += "Reference Position:   " + std::to_string(startPosition) + eol;
		reportString += "Graph Mapping Offset: " + std::to_string(this->m_graph_mapping_ptr->position) + eol;
		reportString += "---------------------------------------------------------" + eol;

		return reportString;
	}
}
}
