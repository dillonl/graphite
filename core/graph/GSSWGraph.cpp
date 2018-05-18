#include "GSSWGraph.h"
#include "core/alignment/AlignmentReporter.h"

#include <mutex>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <thread>

namespace graphite
{

	uint32_t GSSWGraph::s_next_id = 0;
	std::mutex GSSWGraph::s_lock;
	GSSWGraph::GSSWGraph(IReference::SharedPtr referencePtr, VariantList::SharedPtr variantListPtr, Region::SharedPtr regionPtr, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue) :
		m_reference_ptr(referencePtr),
		m_match(matchValue),
		m_mismatch(misMatchValue),
		m_gap_open(gapOpenValue),
		m_gap_extension(gapExtensionValue),
		m_variant_list_ptr(variantListPtr),
		m_region_ptr(regionPtr),
		m_total_graph_length(0),
		m_skipped(false)
	{
		gssw_sse2_disable();
		this->m_nt_table = gssw_create_nt_table();
		this->m_mat = gssw_create_score_matrix(this->m_match, this->m_mismatch);
		this->m_graph_ptr = gssw_graph_create(30);
	}

	GSSWGraph::~GSSWGraph()
	{
		gssw_graph_destroy(this->m_graph_ptr);
		free(this->m_nt_table);
		free(this->m_mat);
	}

	void GSSWGraph::constructGraph()
	{
		int64_t referenceSize;
		std::vector< gssw_node* > altAndRefVertices;
		position currentReferencePosition = this->m_region_ptr->getStartPosition();

		std::vector< IVariant::SharedPtr > variantPtrs = this->m_variant_list_ptr->getAllVariantPtrs();
		for (IVariant::SharedPtr variantPtr : variantPtrs)
		{
			if (variantPtr->shouldSkip())
			{
				m_skipped = true;
				continue;
			}
			referenceSize = variantPtr->getPosition() - currentReferencePosition;
			if (referenceSize > 0)
			{
				auto refRegionPtr = std::make_shared< Region >(this->m_region_ptr->getReferenceID(), currentReferencePosition, variantPtr->getPosition() - 1, Region::BASED::ONE); // minus one because we don't want to include the actual variant position
				std::string referenceSequenceString = this->m_reference_ptr->getSequenceFromRegion(refRegionPtr);
				auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
				auto referenceNode = addReferenceVertex(variantPtr->getRegions()[0]->getStartPosition(), referenceAllelePtr, altAndRefVertices);
				altAndRefVertices.clear();
				altAndRefVertices.push_back(referenceNode);
				m_total_graph_length += referenceSize;
				currentReferencePosition += referenceSize;
			}

			altAndRefVertices = addAlternateVertices(altAndRefVertices, variantPtr);
			currentReferencePosition += variantPtr->getReferenceSize();
		}
        referenceSize = this->m_region_ptr->getEndPosition() - currentReferencePosition;
		if (referenceSize > 0)
		{
			auto refRegionPtr = std::make_shared< Region >(this->m_region_ptr->getReferenceID(), currentReferencePosition, this->m_region_ptr->getEndPosition(), Region::BASED::ONE);
			std::string referenceSequenceString = this->m_reference_ptr->getSequenceFromRegion(refRegionPtr);
			auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
			addReferenceVertex(currentReferencePosition, referenceAllelePtr, altAndRefVertices);
		}
	}

	gssw_node* GSSWGraph::addReferenceVertex(position position, IAllele::SharedPtr referenceAllelePtr, std::vector< gssw_node* > altAndRefVertices)
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
		size_t tmpLength = 0;
	    std::vector< gssw_node* > vertices;
		for (auto altAllelePtr : variantPtr->getAltAllelePtrs())
		{
			auto altAlleleNode = gssw_node_create_alt(variantPtr->getPosition(), variantPtr->getRefAllelePtr()->getSequence(), variantPtr->getRefAllelePtr()->getLength(), altAllelePtr, false, this->m_nt_table, this->m_mat);
			gssw_graph_add_node(this->m_graph_ptr, altAlleleNode);
			vertices.emplace_back(altAlleleNode);
			if (altAllelePtr->getLength() > tmpLength) { tmpLength = altAllelePtr->getLength(); }
		}

		gssw_node* variantReferenceNode = gssw_node_create_alt(variantPtr->getPosition(), variantPtr->getRefAllelePtr()->getSequence(), variantPtr->getRefAllelePtr()->getLength(), variantPtr->getRefAllelePtr(), true, this->m_nt_table, this->m_mat);
		if (variantPtr->getRefAllelePtr()->getLength() > tmpLength) { tmpLength = variantPtr->getRefAllelePtr()->getLength(); }

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
		gssw_graph_fill(this->m_graph_ptr, alignmentPtr->getSequence(), this->m_nt_table, this->m_mat, this->m_gap_open, this->m_gap_extension, 0, 0, 15, 2, true);
		gssw_graph_mapping* graphMapping = gssw_graph_trace_back(this->m_graph_ptr,alignmentPtr->getSequence(),alignmentPtr->getLength(),this->m_nt_table,this->m_mat,m_gap_open,m_gap_extension,0,0);

		auto graphMappingDeletor = [](gssw_graph_mapping* gm)
		{
			gssw_graph_mapping_destroy(gm);
		};
		return std::shared_ptr< gssw_graph_mapping >(graphMapping, graphMappingDeletor);
	}

	void GSSWGraph::recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr)
	{
		throw "you will need to implement VariantList::rewind";
		auto alignmentReport = std::make_shared< AlignmentReport >(this->m_reference_ptr, this->m_variant_list_ptr, alignmentPtr, graphMapping, this->m_region_ptr->getStartPosition());
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

	void getAllPaths(gssw_node* node, std::string currentPath, std::string nodeIDs, int numberOfSibs, std::vector< std::tuple< std::string, std::string > >& paths)
	{
		std::string tmpPath = std::string(node->seq, node->len);
		if (numberOfSibs == 0)
		{
			std::transform(tmpPath.begin(), tmpPath.end(), tmpPath.begin(), tolower);
		}
		else
		{
			std::transform(tmpPath.begin(), tmpPath.end(), tmpPath.begin(), ::toupper);
		}
		currentPath += tmpPath;
		nodeIDs += (node->id % 2 == 0) ? ":REF:" : ":ALT:";
		if (node->count_next == 0)
		{
			paths.emplace_back(std::make_tuple(currentPath, nodeIDs));
		}
		else
		{
			for (auto i = 0; i < node->count_next; ++i)
			{
				getAllPaths(node->next[i], currentPath, nodeIDs, node->count_next - 1, paths);
			}
		}
	}

	std::vector< std::tuple< std::string, std::string > > GSSWGraph::generateAllPaths()
	{
		std::vector< std::tuple< std::string, std::string > > paths;
		getAllPaths(this->m_graph_ptr->nodes[0], "", "", 0, paths);
		return paths;
	}

}
