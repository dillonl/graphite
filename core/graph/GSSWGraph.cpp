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
	GSSWGraph::GSSWGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, Region::SharedPtr regionPtr, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue, uint32_t numGraphCopies) :
		IGraph(referencePtr, variantListPtr),
		m_match(matchValue),
		m_mismatch(misMatchValue),
		m_gap_open(gapOpenValue),
		m_gap_extension(gapExtensionValue),
		m_variant_list_ptr(variantListPtr),
		m_region_ptr(regionPtr),
		m_total_graph_length(0),
		m_skipped(false),
		m_num_graph_copies(numGraphCopies)
	{
		this->m_nt_table = gssw_create_nt_table();
		this->m_mat = gssw_create_score_matrix(this->m_match, this->m_mismatch);
		this->m_graph_ptr = gssw_graph_create(100);
	}

	GSSWGraph::~GSSWGraph()
	{
	}

	void GSSWGraph::constructGraph()
	{
        //std::cout << "Entered constructGraph" << std::endl;
		int64_t referenceSize;
		IVariant::SharedPtr variantPtr = nullptr;
		std::vector< gssw_node* > altAndRefVertices;
		position currentReferencePosition = this->m_region_ptr->getStartPosition();
        //std::cout << "currentReferencePosition: " << currentReferencePosition << std::endl;

		while (this->m_variant_list_ptr->getNextVariant(variantPtr))
		{
            //std::cout << "Entered while loop." << std::endl;
            // Instantiate vectors for the sequence headers and related sequences.
            std::vector< std::string > sequenceHeaders(variantPtr->getAltAllelePtrs().size() + 1);
            std::vector< std::string > sequences(variantPtr->getAltAllelePtrs().size() + 1);

			if (variantPtr->shouldSkip())
			{
				m_skipped = true;
				continue;
			}
			referenceSize = variantPtr->getPosition() - currentReferencePosition;
            //std::cout << "referenceSize: " << referenceSize << std::endl;

            // refSeqEnd should be calculated correctly if the alt sequence region is the same length as the ref sequences.
            std::string refSeqStart = std::string(currentReferencePosition + this->m_reference_ptr->getSequence(), referenceSize);
            std::string refSeqEnd = std::string(currentReferencePosition + this->m_reference_ptr->getSequence() + variantPtr->getRefAllelePtr()->getLength(), referenceSize);

			if (referenceSize > 0)
			{
                //std::cout << "Entered referenceSize if statement." << std::endl;
				auto refRegionPtr = std::make_shared< Region >(this->m_region_ptr->getReferenceID(), currentReferencePosition, variantPtr->getPosition() - 1, Region::BASED::ONE); // minus one because we don't want to include the actual variant position
                //std::cout << "currentReferencePosition: " << currentReferencePosition << std::endl;
                //std::cout << "Setup refRegionPtr." << std::endl;
                //std::cout << "getRegionString(): " << refRegionPtr->getRegionString() << std::endl;
				std::string referenceSequenceString = this->m_reference_ptr->getSequenceFromRegion(refRegionPtr);
                //std::cout << "Setup referenceSequenceString." << std::endl;
				auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
                // std::cout << "Setup referenceAllelePtr." << std::endl;
				auto referenceNode = addReferenceVertex(variantPtr->getRegions()[0]->getStartPosition(), referenceAllelePtr, altAndRefVertices);
                // std::cout << "Setup referenceNode. " << std::endl;
				altAndRefVertices.clear();
				altAndRefVertices.push_back(referenceNode);
                // std::cout << "Added referenceNode." << std::endl;
				m_total_graph_length += referenceSize;
				currentReferencePosition += referenceSize;
			}
            //std::cout << "Exited referenceSize if statement. " << std::endl;

            // Construct paths by Variant.
            // Create sequence headers and insert the starting reference sequence.
            for (int i = 0; i < sequenceHeaders.size(); i++)
            {
                sequenceHeaders[i] = ">" + variantPtr->getChrom() + ":" + std::to_string(variantPtr->getPosition()) + ":" + std::to_string(i);
                sequences[i] = refSeqStart;
            }

            // Append the rest of the reference sequence to refSeqStart.
            sequences[0] += variantPtr->getRefAllelePtr()->getSequence() + refSeqEnd;

            // Add the alternate alleles as well as refSeqEnd. Need to make sure this is done correctly.
            for (int i = 0; i < variantPtr->getAltAllelePtrs().size(); ++i)
            {
                sequences[i + 1] += variantPtr->getAltAllelePtrs()[i]->getSequence() + refSeqEnd;
            }

            // Append headers and paths to member variable.
            for (int i = 0; i < sequences.size(); ++i)
            {
                m_graph_path_headers.emplace_back(sequenceHeaders[i]);
                m_graph_path_sequences.emplace_back(sequences[i]);
            }

			altAndRefVertices = addAlternateVertices(altAndRefVertices, variantPtr);
			currentReferencePosition += variantPtr->getReferenceSize();
		}
        //std::cout << "Exited while loop." << std::endl;
        referenceSize = this->m_region_ptr->getEndPosition() - currentReferencePosition;
		if (referenceSize > 0)
		{
			auto refRegionPtr = std::make_shared< Region >(this->m_region_ptr->getReferenceID(), currentReferencePosition, this->m_region_ptr->getEndPosition(), Region::BASED::ONE);
			std::string referenceSequenceString = this->m_reference_ptr->getSequenceFromRegion(refRegionPtr);
			auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
			addReferenceVertex(currentReferencePosition, referenceAllelePtr, altAndRefVertices);
		}
		generateGraphCopies();
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

	GSSWGraph::GSSWGraphMappingPtr GSSWGraph::traceBackAlignment(IAlignment::SharedPtr alignmentPtr, std::shared_ptr< GSSWGraphContainer > graphContainer)
	{
		gssw_graph* g = graphContainer->graph_ptr;
		int8_t* nt_table = graphContainer->nt_table;
		int8_t* mat = graphContainer->mat;

		gssw_graph_fill(g, alignmentPtr->getSequence(), alignmentPtr->getLength(), nt_table, mat, this->m_gap_open, this->m_gap_extension, 15, 2);
		gssw_graph_mapping* graphMapping = gssw_graph_trace_back(g, alignmentPtr->getSequence(), alignmentPtr->getLength(),m_match,m_mismatch,m_gap_open,m_gap_extension);

		{
			std::unique_lock< std::mutex > lock(m_traceback_lock);
			m_graph_container_ptrs_queue.emplace(graphContainer);
		}
		this->m_condition.notify_one();

		gssw_node_cigar* nc = graphMapping->cigar.elements;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			nc->node->cigar = nc->cigar;
		}

		auto graphMappingDeletor = [](gssw_graph_mapping* gm)
		{
			gssw_graph_mapping_destroy(gm);
		};
		return std::shared_ptr< gssw_graph_mapping >(graphMapping, graphMappingDeletor);
	}

	void GSSWGraph::recordAlignmentVariants(std::shared_ptr< gssw_graph_mapping > graphMapping, IAlignment::SharedPtr alignmentPtr)
	{
		 // this->m_variant_list_ptr->rewind();
		throw "you will need to implement IVariantList::rewind";
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

	std::shared_ptr< GSSWGraphContainer > GSSWGraph::getGraphContainer()
	{
		std::unique_lock< std::mutex > lock(m_traceback_lock);
		this->m_condition.wait(lock, [this]{ return !this->m_graph_container_ptrs_queue.empty(); });
		auto graphContainerPtr = this->m_graph_container_ptrs_queue.front();
		this->m_graph_container_ptrs_queue.pop();
		return graphContainerPtr;
	}

	void GSSWGraph::generateGraphCopies()
	{
		for (uint32_t tc = 0; tc < m_num_graph_copies; ++tc)
		{
			int8_t* nt_table = gssw_create_nt_table();
			int8_t* mat = mat = gssw_create_score_matrix(this->m_match, this->m_mismatch);
			gssw_graph* g = gssw_graph_create(100);

			std::unordered_map< int, gssw_node* > oldToNewNodeMap;
			for (auto i = 0; i < m_graph_ptr->size; ++i)
			{
				auto node = gssw_node_copy(this->m_graph_ptr->nodes[i], nt_table);
				gssw_graph_add_node(g, node);

				oldToNewNodeMap.emplace(this->m_graph_ptr->nodes[i]->id, node);
			}
			for (auto i = 0; i < m_graph_ptr->size; ++i)
			{
				gssw_node* oldStartNode = this->m_graph_ptr->nodes[i];
				for (auto nextIdx = 0; nextIdx < this->m_graph_ptr->nodes[i]->count_next; ++nextIdx)
				{
					gssw_node* oldEndNode = this->m_graph_ptr->nodes[i]->next[nextIdx];
					if (oldToNewNodeMap.find(oldStartNode->id) == oldToNewNodeMap.end() || oldToNewNodeMap.find(oldEndNode->id) == oldToNewNodeMap.end()) { std::cout << "skipping!!!" << std::endl; }
					gssw_node* newStartNode = oldToNewNodeMap[oldStartNode->id];
					gssw_node* newEndNode = oldToNewNodeMap[oldEndNode->id];
					gssw_nodes_add_edge(newStartNode, newEndNode);
				}
			}
			auto graphContainerPtr = std::make_shared< GSSWGraphContainer >(nt_table, mat, g);
			m_graph_container_ptrs.emplace_back(graphContainerPtr);
			m_graph_container_ptrs_queue.emplace(graphContainerPtr);
		}
		auto graphContainerPtr = std::make_shared< GSSWGraphContainer >(this->m_nt_table, this->m_mat, this->m_graph_ptr);
		m_graph_container_ptrs.emplace_back(graphContainerPtr);
		m_graph_container_ptrs_queue.emplace(graphContainerPtr);
	}

    std::vector< std::string > GSSWGraph::getGraphPathHeaders()
    {
        return m_graph_path_headers;
    }

    std::vector< std::string > GSSWGraph::getGraphPathSequences()
    {
        return m_graph_path_sequences;
    }
}
