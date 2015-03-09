#include "GSSWGraph.h"
#include "vg/graph/ReferenceNode.h"
#include "vg/graph/SNPNode.h"
#include "vg/graph/IVariantNode.h"

#include <iostream>
#include <fstream>
#include <map>

namespace gwiz
{
namespace gssw
{
	bool GSSWGraph::PrintStuff = false;

	GSSWGraph::GSSWGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignmentReader::SharedPtr alignmentReader) :
		IGraph(referencePtr, variantListPtr), m_alignment_reader(alignmentReader), m_match(2), m_mismatch(2), m_gap_open(3), m_gap_extension(1)

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
		position alignmentReaderStartPosition = this->m_alignment_reader->getRegion()->getStartPosition();
		position alignmentReaderEndPosition = this->m_alignment_reader->getRegion()->getEndPosition();
		// size_t readLength = this->m_alignment_reader->getAverageReadLength();
		position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
		size_t referenceOffset = (startPosition < alignmentReaderStartPosition) ? (alignmentReaderStartPosition - startPosition) : 0;
		int64_t referenceSize;
		Variant::SharedPtr variantPtr;
		std::vector< gssw_node* > altAndRefVertices;
		size_t graphSize = 0;
		int count = 0;
		static std::mutex mutex;
		while (getNextCompoundVariant(variantPtr))
		{
			referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
			if (referenceSize > 0)
			{
				auto referenceNode = gssw_node_create_alt(this->m_reference_ptr->getSequence() + referenceOffset, referenceSize, this->m_nt_table, this->m_mat);
			    addReference(altAndRefVertices, referenceNode);
				altAndRefVertices.clear();
				altAndRefVertices.push_back(referenceNode);
			}
			size_t variantReferenceSize;
			altAndRefVertices = addAlternateVertices(altAndRefVertices, variantPtr, variantReferenceSize);
			referenceOffset += referenceSize + variantReferenceSize;
		}
		position endPosition = (this->m_reference_ptr->getRegion()->getEndPosition() > alignmentReaderEndPosition) ? alignmentReaderEndPosition : this->m_reference_ptr->getRegion()->getEndPosition();
		referenceSize = endPosition - (startPosition + referenceOffset);
		if (referenceSize > 0)
		{
			auto referenceNode = gssw_node_create_alt(this->m_reference_ptr->getSequence() + referenceOffset, referenceSize, this->m_nt_table, this->m_mat);
			addReference(altAndRefVertices, referenceNode);
		}
		graphConstructed();
	}

	gssw_node* GSSWGraph::addReference(std::vector< gssw_node* > altAndRefVertices, gssw_node* referenceNode)
	{
		gssw_graph_add_node(this->m_graph_ptr, referenceNode);
		for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
		{
			gssw_nodes_add_edge((*iter), referenceNode);
		}
		return referenceNode;
	}

	std::vector< gssw_node* > GSSWGraph::addAlternateVertices(std::vector< gssw_node* > altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize)
	{
	    std::vector< gssw_node* > vertices;
		for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
		{
			INode::SharedPtr variantNodePtr = vg::IVariantNode::BuildVariantNodes(variantPtr, i);
			vertices.push_back(addGSSWAlternateNode(variantNodePtr));
		}
		size_t referenceOffset = variantPtr->getPosition() - this->m_reference_ptr->getRegion()->getStartPosition();
		gssw_node* variantReferenceNode = gssw_node_create_alt(this->m_reference_ptr->getSequence() + referenceOffset, variantPtr->getRef()[0].size(), this->m_nt_table, this->m_mat);
		gssw_graph_add_node(this->m_graph_ptr, variantReferenceNode);
		vertices.push_back(variantReferenceNode);
		variantReferenceSize = variantPtr->getRef()[0].size();
		for (auto parentVertexIter = altAndRefVertices.begin(); parentVertexIter != altAndRefVertices.end(); ++parentVertexIter)
		{
			for (auto childVertexIter = vertices.begin(); childVertexIter != vertices.end(); ++childVertexIter)
			{
				gssw_nodes_add_edge((*parentVertexIter), (*childVertexIter));
			}
		}
		return vertices;
	}

	std::shared_ptr< gssw_graph_mapping > GSSWGraph::traceBackAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		std::string readSeq = std::string(alignmentPtr->getSequence(), alignmentPtr->getLength());
		gssw_graph_fill(this->m_graph_ptr, readSeq.c_str(), this->m_nt_table, this->m_mat, this->m_gap_open, this->m_gap_extension, 15, 2);
		gssw_graph_mapping* graphMapping = gssw_graph_trace_back (this->m_graph_ptr,readSeq.c_str(),readSeq.size(),m_match,m_mismatch,m_gap_open,m_gap_extension);
		auto graphMappingDeletor = [](gssw_graph_mapping* gm) { gssw_graph_mapping_destroy(gm); };
		return std::shared_ptr< gssw_graph_mapping >(graphMapping, graphMappingDeletor);
	}

	void GSSWGraph::recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr)
	{
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			if (m_node_map.find(nc->node->id) != m_node_map.end())
			{
				auto variantNode = m_node_map[nc->node->id];
				uint32_t counter = 1;
				std::vector< IAlignment::SharedPtr > alignments;
				if (m_variant_counter.find(nc->node->id) != m_variant_counter.end())
				{
					counter = std::get< 1 >(m_variant_counter[nc->node->id]) + 1;
					alignments = std::get< 2 >(m_variant_counter[nc->node->id]);
				}
				alignments.push_back(alignmentPtr);
				m_variant_counter[nc->node->id] = std::make_tuple(variantNode, counter, alignments);
			}
		}
	}

	static std::ofstream vcfCount;

	void GSSWGraph::graphConstructed()
	{
		IAlignment::SharedPtr alignmentPtr;
		while (this->m_alignment_reader->getNextAlignment(alignmentPtr))
		{
			std::shared_ptr< gssw_graph_mapping > graphMapping = traceBackAlignment(alignmentPtr);
			recordAlignmentVariants(graphMapping, alignmentPtr);
		}

		static bool init = false;
		if (!init)
		{
			vcfCount.open("test.txt");
			init = true;
		}
		for (const auto& value : m_variant_counter)
		{
			auto variant = std::get< 0 >(value.second);
			vcfCount << "Variant Count: " << std::get< 1 >(value.second) << " Variant Seq: " << std::string(variant->getSequence(), variant->getLength()) << " " << variant->getPosition() << std::endl;
			for (const auto& alignmentPtr : std::get< 2 >(value.second))
			{
				vcfCount << "Alignment: " << alignmentPtr->getPosition() << std::endl;
			}
		}

	}
}
}
