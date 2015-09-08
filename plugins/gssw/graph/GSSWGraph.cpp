#include "GSSWGraph.h"
#include "vg/graph/ReferenceNode.h"
#include "vg/graph/SNPNode.h"
#include "vg/graph/IVariantNode.h"
#include "core/genotyper/IGenotyper.h"
#include "AlignmentReporter.h"

#include <mutex>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

namespace graphite
{
namespace gssw
{

	GSSWGraph::GSSWGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, position startPosition, size_t graphSize, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue) :
		IGraph(referencePtr, variantListPtr),
		m_match(matchValue),
		m_mismatch(misMatchValue),
		m_gap_open(gapOpenValue),
		m_gap_extension(gapExtensionValue),
		m_variant_list_ptr(variantListPtr),
		m_start_position(startPosition),
		m_graph_size(graphSize)

	{
		this->m_nt_table = gssw_create_nt_table();
		this->m_mat = gssw_create_score_matrix(this->m_match, this->m_mismatch);
		this->m_graph_ptr = gssw_graph_create(100);
	}

	GSSWGraph::~GSSWGraph()
	{
		gssw_graph_destroy(this->m_graph_ptr);
		free(this->m_nt_table);
		free(this->m_mat);
	}

	void GSSWGraph::constructGraph()
	{
		position startPosition = this->m_start_position;
		position endPosition = this->m_start_position + this->m_graph_size;
		size_t graphStartOffset = startPosition - this->m_reference_ptr->getRegion()->getStartPosition();

		size_t referenceOffset = 0;//startPosition;
		int64_t referenceSize;
		IVariant::SharedPtr variantPtr;
		std::vector< gssw_node* > altAndRefVertices;

		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
			// if (!variantPtr->processSV(this->m_reference_ptr)) { continue; }
			auto genotyperVariantPtr = IGenotyper::Instance()->generateVariant(variantPtr->getPosition());
			referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
			if (referenceSize > 0)
			{
				std::string referenceSequenceString = std::string(this->m_reference_ptr->getSequence() + graphStartOffset + referenceOffset, referenceSize);
				auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
				auto referenceNode = addReference((startPosition + referenceOffset), referenceAllelePtr, altAndRefVertices);
				altAndRefVertices.clear();
				altAndRefVertices.push_back(referenceNode);
			}

			altAndRefVertices = addAlternateVertices(altAndRefVertices, variantPtr);
			referenceOffset += referenceSize + variantPtr->getRefAllelePtr()->getLength();
		}

		position currentEndPosition = (this->m_reference_ptr->getRegion()->getEndPosition() > endPosition) ? endPosition : this->m_reference_ptr->getRegion()->getEndPosition();
		referenceSize = currentEndPosition - (startPosition + referenceOffset);
		if (referenceSize > 0)
		{
			std::string referenceSequenceString = std::string(this->m_reference_ptr->getSequence() + graphStartOffset + referenceOffset, referenceSize);
			auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
			addReference((startPosition + referenceOffset), referenceAllelePtr, altAndRefVertices);
		}
	}

	gssw_node* GSSWGraph::addReference(position position, IAllele::SharedPtr referenceAllelePtr, std::vector< gssw_node* > altAndRefVertices)
	{
		this->m_reference_fragments.emplace_back(referenceAllelePtr);
		auto referenceNodePtr = gssw_node_create_alt(position, referenceAllelePtr->getSequence(), referenceAllelePtr->getLength(), referenceAllelePtr, true, this->m_nt_table, this->m_mat);
		gssw_graph_add_node(this->m_graph_ptr, referenceNodePtr);
		for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
		{
			gssw_nodes_add_edge((*iter), referenceNodePtr);
		}
		return referenceNodePtr;
	}

	std::vector< gssw_node* > GSSWGraph::addAlternateVertices(const std::vector< gssw_node* >& altAndRefVertices, IVariant::SharedPtr variantPtr)
	{
		size_t referenceOffset = variantPtr->getPosition() - this->m_reference_ptr->getRegion()->getStartPosition();
	    std::vector< gssw_node* > vertices;
		for (auto altAllelePtr : variantPtr->getAltAllelePtrs())
		{

			auto altAlleleNode = gssw_node_create_alt(variantPtr->getPosition(), variantPtr->getRefAllelePtr()->getSequence(), variantPtr->getRefAllelePtr()->getSequencePtr()->getLength(), altAllelePtr, false, this->m_nt_table, this->m_mat);
			gssw_graph_add_node(this->m_graph_ptr, altAlleleNode);
			vertices.emplace_back(altAlleleNode);
		}

		gssw_node* variantReferenceNode = gssw_node_create_alt(variantPtr->getPosition(), variantPtr->getRefAllelePtr()->getSequence(), variantPtr->getRefAllelePtr()->getLength(), variantPtr->getRefAllelePtr(), true, this->m_nt_table, this->m_mat);

		gssw_graph_add_node(this->m_graph_ptr, variantReferenceNode);
		vertices.push_back(variantReferenceNode);
		for (auto parentNode : altAndRefVertices)
		{
			for (auto childNode : vertices)
			{
				gssw_nodes_add_edge(parentNode, childNode);
			}
		}
		return vertices;
	}

	GSSWGraph::GSSWGraphMappingPtr GSSWGraph::traceBackAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		gssw_graph_fill(this->m_graph_ptr, alignmentPtr->getSequence(), alignmentPtr->getLength(), this->m_nt_table, this->m_mat, this->m_gap_open, this->m_gap_extension, 15, 2);
		gssw_graph_mapping* graphMapping = gssw_graph_trace_back(this->m_graph_ptr, alignmentPtr->getSequence(), alignmentPtr->getLength(),m_match,m_mismatch,m_gap_open,m_gap_extension);
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			nc->node->cigar = nc->cigar;
		}
		// gssw_print_graph_mapping(graphMapping);
		// gssw_graph_print_score_matrices(this->m_graph_ptr, alignmentPtr->getSequence(), alignmentPtr->getLength());
		auto graphMappingDeletor = [](gssw_graph_mapping* gm) { gssw_graph_mapping_destroy(gm); };
		return std::shared_ptr< gssw_graph_mapping >(graphMapping, graphMappingDeletor);
	}

	void GSSWGraph::recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr)
	{
		 // this->m_variant_list_ptr->rewind();
		throw "you will need to implement IVariantList::rewind";
		 auto alignmentReport = std::make_shared< AlignmentReport >(this->m_reference_ptr, this->m_variant_list_ptr, alignmentPtr, graphMapping, this->m_start_position);
		 AlignmentReporter::Instance()->addAlignmentReport(alignmentReport);
	}

	IVariant::SharedPtr GSSWGraph::getVariantFromNodeID(const uint32_t nodeID)
	{
		if (this->m_variants_map.find(nodeID) != this->m_variants_map.end())
		{
			return this->m_variants_map[nodeID];
		}
		return nullptr;
	}

	void GSSWGraph::graphConstructed()
	{
	}

	IAllele::SharedPtr GSSWGraph::getAllelePtrFromNodeID(uint32_t id)
	{
		auto iter = m_node_id_to_allele_ptrs.find(id);
		if (iter != m_node_id_to_allele_ptrs.end())
		{
			return iter->second;
		}
		else
		{
			return nullptr;
		}
	}
}
}
