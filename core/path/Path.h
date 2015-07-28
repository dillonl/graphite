#ifndef GWIZ_PATH_H
#define GWIZ_PATH_H

#include "IPath.h"

#include <vector>
#include <memory>

namespace gwiz
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

	private:
		void computeAndSetHash();

		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		IAlignment::SharedPtr m_alignment_ptr;
		uint32_t m_sw_percentage;
		size_t m_hash;
	};
}

#endif //GWIZ_PATH_H
