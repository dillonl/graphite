#ifndef GWIZ_PATHMANAGER_H
#define GWIZ_PATHMANAGER_H

#include "core/alignment/IAlignment.h"
#include "core/allele/IAllele.h"
#include "core/path/IPath.h"

#include <boost/noncopyable.hpp>

#include <unordered_map>
#include <vector>
#include <mutex>

namespace gwiz
{
	class PathManager : private boost::noncopyable
	{
	public:
		static PathManager* Instance();
		IPath::SharedPtr getPath(size_t hash);
		std::vector< IPath::SharedPtr > getPaths(IAllele::SharedPtr allelePtr);
		std::vector< IPath::SharedPtr > getPaths(IAlignment::SharedPtr alignmentPtr);

		void addPath(IAlignment::SharedPtr, IPath::SharedPtr);

	private:
		PathManager() {}
		~PathManager() {}

		std::unordered_map< size_t, IPath::SharedPtr > m_hash_path;
		std::unordered_map< IAllele::SharedPtr, std::vector< IPath::SharedPtr > > m_allele_path_ptrs;
		std::unordered_map< IAlignment::SharedPtr, std::vector< IPath::SharedPtr > > m_alignment_path_ptrs;
		std::mutex m_path_mutex;

		static PathManager* s_instance;
	};
}

#endif //GWIZ_PATHMANAGER_H
