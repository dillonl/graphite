/*
#ifndef GRAPHITE_VISUALIZATIONTOOLKIT_H
#define GRAPHITE_VISUALIZATIONTOOLKIT_H

#include "core/alignment/IAlignment.h"
#include "core/mapping/GSSWMapping.h"
#include "core/file/FastaFileWriter.h"
#include "core/file/BamFileWriter.h"
#include "core/graph/GSSWGraph.h"
#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class VisualizationToolKit : private Noncopyable
	{
	public:
		typedef std::shared_ptr< VisualizationToolKit > SharedPtr;
		VisualizationToolKit(const std::string& outputPath, const std::vector< std::string >& bamPaths, int matchValue);
		virtual ~VisualizationToolKit();

		void setAlignmentAndMapping(IAlignment::SharedPtr alignmentPtr, GSSWGraph::SharedPtr gsswGraphPtr, GSSWMapping::SharedPtr refMapping, GSSWMapping::SharedPtr altMapping);
		void closeResources();

	private:
		void initializeKit(const std::vector< std::string >& bamPaths);
		std::vector< BamTools::CigarOp > cigarToBamToolsCigar(const std::vector< std::tuple< char, uint32_t > >& cigar);

		int m_match_value;
		std::string m_output_path;
		std::vector< BamFileWriter::SharedPtr > m_bam_writer_ptrs;
		std::vector< FastaFileWriter::SharedPtr > m_fasta_writer_ptrs;
		std::unordered_set< std::string > m_descriptions;
		std::unordered_map< std::string, FastaFileWriter::SharedPtr > m_fasta_files;
		std::unordered_map< std::string, BamFileWriter::SharedPtr > m_bam_files;
	};
}

#endif // GRAPHITE_VISUALIZATIONTOOLKIT_H
*/
