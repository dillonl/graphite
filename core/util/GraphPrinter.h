#ifndef GRAPHITE_GRAPHPRINTER_H
#define GRAPHITE_GRAPHPRINTER_H

#include "core/util/Noncopyable.hpp"
#include "core/graph/Node.h"
#include "core/alignment/Alignment.h"

#include "gssw.h"

#include <tuple>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace graphite
{
	class Graph;
	class GraphPrinter : private Noncopyable
	{
	public:
		typedef std::shared_ptr< GraphPrinter > SharedPtr;
		GraphPrinter(Graph* graphPtr);
		~GraphPrinter();

		void registerTraceback(gssw_graph_mapping* graphMapping, Alignment::SharedPtr alignmentPtr, float sswScore);
		void registerUnalignedRead(Alignment::SharedPtr alignmentPtr, std::string graphCigarString, float sswScore);
		void printGraph();

	private:
		struct MappingContainer
		{
			uint32_t m_ssw_score;
			std::string m_sequence;
			std::string m_read_name;
			std::string m_cigar_str;
			std::string m_graph_cigar;
			/* std::vector< std::tuple< uint32_t, char > > m_cigar; */
			uint32_t m_position_offset;
		};
		struct UnMappedContainer
		{
			uint32_t m_ssw_score;
			std::string m_sequence;
			std::string m_read_name;
			std::string m_cigar_str;
			std::string m_graph_cigar;
			uint32_t m_position;
		};
		std::shared_ptr< std::unordered_set< uint32_t > > getKey(std::unordered_set< uint32_t >& set);
		uint32_t getGraphOffset(Node* nodePtr);

		/* Graph* m_graph_ptr; */
		std::unordered_map< std::shared_ptr< std::unordered_set< uint32_t > >, std::vector< Node::SharedPtr > > m_graph_path_node_ids;
		std::unordered_map< std::shared_ptr< std::unordered_set< uint32_t > >, std::vector< std::shared_ptr< MappingContainer > > > m_graph_mapping_container_ptrs;
		std::unordered_set< std::string > m_alignment_readname_tracker_map;
		std::vector< std::shared_ptr< UnMappedContainer > > m_unmapped_read_containers;
	};
}

#endif //GRAPHITE_GRAPHPRINTER_H
