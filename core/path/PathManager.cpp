#include "PathManager.h"

namespace gwiz
{
	PathManager* PathManager::s_instance = nullptr;

	PathManager* PathManager::Instance()
	{
		if (s_instance == nullptr)
		{
			s_instance = new PathManager();
		}
		return s_instance;
	}

	IPath::SharedPtr PathManager::getPath(size_t hash)
	{
		std::lock_guard< std::mutex > lock(this->m_path_mutex);
		auto iter = this->m_hash_path.find(hash);
		if (iter != this->m_hash_path.end())
		{
			return iter->second;
		}
		else
		{
			return nullptr;
		}
	}

	std::vector< IPath::SharedPtr > PathManager::getPaths(IAllele::SharedPtr allelePtr)
	{
		std::lock_guard< std::mutex > lock(this->m_path_mutex);
		auto iter = this->m_allele_path_ptrs.find(allelePtr);
		if (iter != this->m_allele_path_ptrs.end())
		{
			return iter->second;
		}
		else
		{
			std::vector< IPath::SharedPtr > pathPtrs;
			return pathPtrs;
		}
	}

	std::vector< IPath::SharedPtr > PathManager::getPaths(IAlignment::SharedPtr alignmentPtr)
	{
		std::lock_guard< std::mutex > lock(this->m_path_mutex);
		auto iter = this->m_alignment_path_ptrs.find(alignmentPtr);
		if (iter != this->m_alignment_path_ptrs.end())
		{
			return iter->second;
		}
		else
		{
			std::vector< IPath::SharedPtr > pathPtrs;
			return pathPtrs;
		}
	}

	void PathManager::addPath(IAlignment::SharedPtr alignmentPtr, IPath::SharedPtr pathPtr)
	{
		std::lock_guard< std::mutex > lock(this->m_path_mutex);
		auto iter = this->m_alignment_path_ptrs.find(alignmentPtr);
		if (iter != this->m_alignment_path_ptrs.end())
		{
			for (auto pathPtr : iter->second)
			{

			}
		}
	}
}
