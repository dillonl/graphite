#include "GraphPrinter.h"

#include "core2/graph/Graph.h"

#include <algorithm>

namespace graphite
{
	GraphPrinter::GraphPrinter(Graph* graphPtr) :
		m_graph_ptr(graphPtr)
	{
		auto paths = graphPtr->generateAllPaths();
		for (auto path : paths)
		{
			std::shared_ptr< std::unordered_set< uint32_t > > pathNodeIDs = std::make_shared< std::unordered_set< uint32_t > >();
			for (auto nodePtr : path)
			{
				pathNodeIDs->emplace(nodePtr->getID());
			}
			m_graph_path_node_ids.emplace(pathNodeIDs, path);
		}
	}

	GraphPrinter::~GraphPrinter()
	{
	}

	void GraphPrinter::registerTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, float sswScore)
	{
		static std::mutex l;
		l.lock();
		std::shared_ptr< MappingContainer > mappingContainerPtr = std::make_shared< MappingContainer >();
		std::vector< std::string > pathNodeIDs;
		std::vector< std::tuple< uint32_t, char > > cigar;
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		char prevType = 0;
		uint32_t prevLength = 0;
		std::unordered_set< uint32_t > nodeIDs;
		Node* firstNodePtr = nullptr;
		uint32_t softclipOffset = 0;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			gssw_node* gsswNode = graphMapping->cigar.elements[i].node;
			uint32_t nodeID = gsswNode->id;
			Node* nodePtr = (Node*)gsswNode->data;
			if (firstNodePtr == nullptr)
			{
				firstNodePtr = nodePtr;
			}
			nodeIDs.emplace(nodeID);
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				if (i == 0 && j == 0 && nc->cigar->elements[0].type == 'S')
				{
					softclipOffset = nc->cigar->elements[0].length;
				}
				if (j == 0 && prevType == nc->cigar->elements[0].type)
				{
					mappingContainerPtr->m_cigar_str += std::to_string(nc->cigar->elements[0].length + prevLength) + std::string(1, nc->cigar->elements[0].type);
				}
				else if (j != nc->cigar->length - 1)
				{
					mappingContainerPtr->m_cigar_str += std::to_string(nc->cigar->elements[0].length) + std::string(1, nc->cigar->elements[0].type);
				}
				mappingContainerPtr->m_graph_cigar += std::to_string(nc->cigar->elements[j].length);
				mappingContainerPtr->m_graph_cigar += std::string(1, nc->cigar->elements[j].type);
				prevType = nc->cigar->elements[j].type;
				prevLength = nc->cigar->elements[j].length;
				// cigar.emplace_back(std::make_tuple< uint32_t, char >(nc->cigar->elements[j].length, nc->cigar->elements[j].type));
			}
		}
		mappingContainerPtr->m_ssw_score = sswScore;
		mappingContainerPtr->m_cigar_str += std::to_string(prevLength) + std::string(1, prevType);
		std::shared_ptr< std::unordered_set< uint32_t > > pathKey = getKey(nodeIDs);
		uint32_t startPosition = getGraphOffset(firstNodePtr);
		mappingContainerPtr->m_position_offset = graphMapping->position + startPosition;
		// std::cout << mappingContainerPtr->m_position_offset << std::endl;
		if (pathKey == nullptr)
		{
			std::cout << "no key" << std::endl;
			return;
		}

		mappingContainerPtr->m_sequence = bamAlignmentPtr->QueryBases;
		mappingContainerPtr->m_read_name = bamAlignmentPtr->Name;
		std::cout << bamAlignmentPtr->Name << std::endl;
		if (softclipOffset >= 0)
		{

			mappingContainerPtr->m_sequence.erase(0, softclipOffset);
		}

		auto iter = m_graph_mapping_container_ptrs.find(pathKey);
		if (iter == m_graph_mapping_container_ptrs.end())
		{
			std::vector< std::shared_ptr< MappingContainer > > mappingContainerPtrs;
			mappingContainerPtrs.emplace_back(mappingContainerPtr);
			m_graph_mapping_container_ptrs.emplace(pathKey, mappingContainerPtrs);
			iter = m_graph_mapping_container_ptrs.find(pathKey);
		}
		else
		{
			std::vector< std::shared_ptr< MappingContainer > >* tmp = &iter->second;
			tmp->emplace_back(mappingContainerPtr);
		}
		l.unlock();
	}

	uint32_t GraphPrinter::getGraphOffset(Node* nodePtr)
	{
		uint32_t offset = 0;
		Node* currentPtr = nodePtr;
		while (currentPtr->getReferenceInNode() != nullptr)
		{
			if (currentPtr->getReferenceInNode() != nullptr)
			{
				currentPtr = currentPtr->getReferenceInNode().get();
				offset += currentPtr->getSequence().size();
			}
			else
			{
				currentPtr = nullptr;
			}
		}
		return offset;
	}

	std::shared_ptr< std::unordered_set< uint32_t > > GraphPrinter::getKey(std::unordered_set< uint32_t >& set)
	{
		for (auto iter : m_graph_path_node_ids)
		{
			int nodeCount = 0;
			for (Node::SharedPtr nodePtr : iter.second)
			{
				if (set.find(nodePtr->getID()) != set.end())
				{
					++nodeCount;
				}
			}
			if (nodeCount == set.size())
			{
				return iter.first;
			}
		}
		return nullptr;
	}

	// TODO:
	// create print all paths function that prints all paths with the metadata of each alignment
	void GraphPrinter::printGraph()
	{
		for (auto iter : m_graph_path_node_ids)
		{
			std::vector< std::string > output;
			std::shared_ptr< std::unordered_set< uint32_t > > nodeIDKey = iter.first;
			std::vector< Node::SharedPtr > nodePath = iter.second;
			auto graphContainerIter = m_graph_mapping_container_ptrs.find(nodeIDKey);
			if (graphContainerIter == m_graph_mapping_container_ptrs.end())
			{
				continue;
			}
			std::string graphPath = "";
			std::string graphNodeID = "";
			for (auto nodePtr : nodePath)
			{
				std::string data = nodePtr->getSequence();
				if (nodePtr->getInNodes().size() == 1 && nodePtr->getOutNodes().size() == 1)
				{
					std::transform(data.begin(), data.end(), data.begin(), ::toupper);
				}
				else
				{
					std::transform(data.begin(), data.end(), data.begin(), ::tolower);
				}
				graphPath += data;
				graphNodeID += (nodePtr->getAlleleType() == Node::ALLELE_TYPE::REF) ? "REF|" : "ALT|";
			}
			std::cout << graphNodeID << std::endl;
			std::cout << graphPath << std::endl;
			std::vector< std::shared_ptr< MappingContainer > > graphMappingContainerPtrs = graphContainerIter->second;
			sort(graphMappingContainerPtrs.begin(), graphMappingContainerPtrs.end(), [](const std::shared_ptr< MappingContainer > lhs, const std::shared_ptr< MappingContainer > rhs) {
					return lhs->m_position_offset < rhs->m_position_offset;
				});

			for (auto graphMappingContainerPtr : graphMappingContainerPtrs)
			{
				std::string line = std::string(graphMappingContainerPtr->m_position_offset, ' ');
				line += graphMappingContainerPtr->m_sequence;
				line += "\t";
				line += std::to_string(graphMappingContainerPtr->m_ssw_score) + "\t" + graphMappingContainerPtr->m_cigar_str + "\t" + graphMappingContainerPtr->m_read_name;

				// output.emplace_back(line);
				std::cout << line << std::endl;

			}
			// for (auto line : output)
			// {
				// std::cout << line << std::endl;
			// }
			// std::cout << std::string(graphPath.size(), '-')<< std::endl;
			// std::cout << std::string(graphPath.size(), '-')<< std::endl;

			std::cout << std::string(400, '-')<< std::endl;
			std::cout << std::string(400, '-')<< std::endl;

		}
	}
}
