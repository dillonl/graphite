#ifndef GRAPHITE_PATH_H
#define GRAPHITE_PATH_H

#include "IPath.h"

#include "gssw.h"

#include <vector>
#include <memory>

namespace graphite
{
class Path : IPath
	{
	public:
		typedef std::shared_ptr< Path > SharedPtr;
        Path();
        ~Path();

		std::vector< IAllele::SharedPtr > getAllelePath() override;
		uint32_t getPathSWPercent() override;
		IAlignment::SharedPtr getAlignment() override;
		size_t getHash() override;

		void addAlleleToPath(IAllele::SharedPtr allelePtr) override;
		void setPathSWPercent(uint32_t swPercent) override;
		void setAlignment(IAlignment::SharedPtr allelePtr) override;
		void setGSSWGraphMapping(std::shared_ptr< gssw_graph_mapping > graphMappingPtr);

		void printLongFormat();

	private:
		void computeAndSetHash();

		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		IAlignment::SharedPtr m_alignment_ptr;
		std::shared_ptr< gssw_graph_mapping > m_graph_mapping_ptr;
		uint32_t m_sw_percentage;
		size_t m_hash;
		std::hash< IAllele* > m_hasher;
	};
}

#endif //GRAPHITE_PATH_H
