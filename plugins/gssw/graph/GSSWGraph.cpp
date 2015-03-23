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
		position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
		size_t referenceOffset = (startPosition < alignmentReaderStartPosition) ? (alignmentReaderStartPosition - startPosition) : 0;
		int64_t referenceSize;
		Variant::SharedPtr variantPtr;
		std::vector< gssw_node* > altAndRefVertices;
		static std::mutex mutex;
		while (getNextCompoundVariant(variantPtr))
		{
			auto genotyperVariantPtr = IGenotyper::Instance()->generateVariant(variantPtr->getPosition());
			referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
			if (referenceSize > 0)
			{
				auto referenceNode = gssw_node_create_alt(this->m_reference_ptr->getSequence() + referenceOffset, referenceSize, this->m_nt_table, this->m_mat);
			    addReference(altAndRefVertices, referenceNode, genotyperVariantPtr);
				altAndRefVertices.clear();
				altAndRefVertices.push_back(referenceNode);
			}
			size_t variantReferenceSize;
			altAndRefVertices = addAlternateVertices(altAndRefVertices, variantPtr, variantReferenceSize, genotyperVariantPtr);
			referenceOffset += referenceSize + variantReferenceSize;
		}
		position endPosition = (this->m_reference_ptr->getRegion()->getEndPosition() > alignmentReaderEndPosition) ? alignmentReaderEndPosition : this->m_reference_ptr->getRegion()->getEndPosition();
		referenceSize = endPosition - (startPosition + referenceOffset);
		if (referenceSize > 0)
		{
			auto referenceNode = gssw_node_create_alt(this->m_reference_ptr->getSequence() + referenceOffset, referenceSize, this->m_nt_table, this->m_mat);
			addReference(altAndRefVertices, referenceNode, nullptr);
		}
		graphConstructed();
	}

	gssw_node* GSSWGraph::addReference(std::vector< gssw_node* > altAndRefVertices, gssw_node* referenceNodePtr, IGenotyperVariant::SharedPtr genotyperVariantPtr)
	{
		gssw_graph_add_node(this->m_graph_ptr, referenceNodePtr);
		for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
		{
			gssw_nodes_add_edge((*iter), referenceNodePtr);
		}
		return referenceNodePtr;
	}

	gssw_node* GSSWGraph::addAlternateNode(INode::SharedPtr variantNodePtr, IGenotyperVariant::SharedPtr genotyperVariantPtr)
	{
		auto variantNode = gssw_node_create_alt(variantNodePtr->getSequence(), variantNodePtr->getLength(), this->m_nt_table, this->m_mat);
		gssw_graph_add_node(this->m_graph_ptr, variantNode);
		auto allelePtr = std::make_shared< GenotyperAllele >(GenotyperAllele::Type::ALTERNATE, variantNode->seq, variantNodePtr->getPosition());
		this->m_genotyper_map[variantNode->id] = allelePtr;
		genotyperVariantPtr->addAllele(allelePtr);
		return variantNode;
	}

	std::vector< gssw_node* > GSSWGraph::addAlternateVertices(std::vector< gssw_node* > altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize, IGenotyperVariant::SharedPtr genotyperVariantPtr)
	{
	    std::vector< gssw_node* > vertices;
		for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
		{
			INode::SharedPtr variantNodePtr = vg::IVariantNode::BuildVariantNodes(variantPtr, i);
			vertices.push_back(addAlternateNode(variantNodePtr, genotyperVariantPtr));
		}
		size_t referenceOffset = variantPtr->getPosition() - this->m_reference_ptr->getRegion()->getStartPosition();
		variantReferenceSize = variantPtr->getRef()[0].size();
		gssw_node* variantReferenceNode = gssw_node_create_alt(this->m_reference_ptr->getSequence() + referenceOffset, variantReferenceSize, this->m_nt_table, this->m_mat);

		auto allelePtr = std::make_shared< GenotyperAllele >(GenotyperAllele::Type::REFERENCE, variantReferenceNode->seq, genotyperVariantPtr->getPosition());
		this->m_genotyper_map[variantReferenceNode->id] = allelePtr;
		genotyperVariantPtr->addAllele(allelePtr);

		gssw_graph_add_node(this->m_graph_ptr, variantReferenceNode);
		vertices.push_back(variantReferenceNode);
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
		gssw_graph_mapping* graphMapping = gssw_graph_trace_back(this->m_graph_ptr,readSeq.c_str(),readSeq.size(),m_match,m_mismatch,m_gap_open,m_gap_extension);
		auto graphMappingDeletor = [](gssw_graph_mapping* gm) { gssw_graph_mapping_destroy(gm); };
		return std::shared_ptr< gssw_graph_mapping >(graphMapping, graphMappingDeletor);
	}

	void GSSWGraph::recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr)
	{
		static std::mutex cigLock;
		cigLock.lock();
		std::string cigarString = "";
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		std::string variantTypeTracebackPositions = "";
		std::string variantTypeTraceback = "";
		std::string genotypeAlleleString = "";
		std::string alignmentMappedString = std::string(alignmentPtr->getSequence(), 100);
		std::string referenceString = std::string((this->m_reference_ptr->getSequence() + alignmentPtr->getPosition()) - (graphMapping->position + 1), 100 + graphMapping->position);
		std::string mappedString = alignmentPtr->isMapped() ? "Mapped" : "Unmapped";
		std::string firstMateString = alignmentPtr->isFirstMate() ? "First" : "Second";
		uint32_t softclipCorrection = 0;
		uint32_t posOffset = 0;
		uint32_t insertionCounter = 0;
		bool firstTimeThrough = true;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				if (i == 0 && j == 0 && nc->cigar->elements[j].type == 'S')
				{
					softclipCorrection = nc->cigar->elements[j].length;
					auto alignmentSize = alignmentMappedString.size();
					if (nc->cigar->elements[j].length == 1)
					{
						softclipCorrection -= 1;
					}
					alignmentMappedString = alignmentMappedString.substr(softclipCorrection, alignmentMappedString.size() - softclipCorrection);
					cigarString += std::to_string(nc->cigar->elements[j].length) + "SC";
				}
				else
				{
					cigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
				}
			}
			std::string suffix = (i < graphMapping->cigar.length - 1) ? "|" : "";
			cigarString += suffix;
			std::string typeString = "R";
			{
				auto genotyperAllele = m_genotyper_map.find(nc->node->id);
				if (genotyperAllele != m_genotyper_map.end())
				{
					if (firstTimeThrough)
					{
						size_t addOne = (i == 0) ? 1 : 0;
						referenceString = std::string((this->m_reference_ptr->getSequence() + genotyperAllele->second->getPosition()) - (genotypeAlleleString.size() + addOne), 100 + graphMapping->position);
						firstTimeThrough = false;
					}
					genotyperAllele->second->addAlignment(alignmentPtr);
					typeString = (genotyperAllele->second->getType() == GenotyperAllele::Type::ALTERNATE) ? "AV" : "RV";
				}
			}
			uint32_t length = (i > 0) ?  nc->node->len : nc->node->len - graphMapping->position;
			if (insertionCounter > 0)
			{
				if (posOffset + (insertionCounter - 1) < alignmentMappedString.size() - 1)
				{
					alignmentMappedString.insert(posOffset + (insertionCounter - 1), " ");
				}
				if (posOffset + graphMapping->position + (insertionCounter - 1) < referenceString.size() - 1)
				{
					referenceString.insert(posOffset + graphMapping->position + (insertionCounter - 1), " ");
				}
			}
			++insertionCounter;
			genotypeAlleleString += std::string(nc->node->seq, nc->node->len) + suffix;
			variantTypeTraceback += typeString + suffix;
			variantTypeTracebackPositions += std::to_string(alignmentPtr->getPosition() + posOffset) + suffix;
			posOffset += length;
		}
		int32_t spaceCount = std::max< int32_t >(graphMapping->position, 0);
		std::string spacing(spaceCount, ' ');
		std::cout << "---------------------------------------------------------" << std::endl;
		std::cout << "Score:              " << graphMapping->score << std::endl;
		std::cout << "Variant Type:       " << variantTypeTraceback << std::endl;
		std::cout << "Variant Positions:  " << variantTypeTracebackPositions << std::endl;
		std::cout << "Cigar String:       " << cigarString << std::endl;
		std::cout << "Reference String:   " << referenceString << std::endl;
		std::cout << "Haplotype String:   " << genotypeAlleleString << std::endl;
		std::cout << "Alignment String:   " << spacing;
		std::cout.write(alignmentMappedString.c_str(), alignmentMappedString.size() - 1);
		std::cout << std::endl;
		std::cout << "Mapping:            " << mappedString << std::endl;
		std::cout << "Mate:               " << firstMateString << std::endl;
		std::cout << "Mapping Position:   " << graphMapping->position << std::endl;
		std::cout << "Alignment Position: " << alignmentPtr->getPosition() << std::endl;
		std::cout << "Reverse Strand:     " << alignmentPtr->isReverseStrand() << std::endl;
		std::cout << "---------------------------------------------------------" << std::endl;
		cigLock.unlock();
	}

	void GSSWGraph::graphConstructed()
	{
		IAlignment::SharedPtr alignmentPtr;
		while (this->m_alignment_reader->getNextAlignment(alignmentPtr))
		{
			if (std::string(alignmentPtr->getSequence()).find("<") != std::string::npos) { continue; } // skip symbolic variants for now
			std::shared_ptr< gssw_graph_mapping > graphMapping = traceBackAlignment(alignmentPtr);
			recordAlignmentVariants(graphMapping, alignmentPtr);
		}
	}
}
}
