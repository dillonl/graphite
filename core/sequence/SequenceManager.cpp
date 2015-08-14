#include "SequenceManager.h"

namespace graphite
{
	SequenceManager* SequenceManager::s_instance = nullptr;

	SequenceManager* SequenceManager::Instance()
	{
		if (s_instance == nullptr)
		{
			s_instance = new SequenceManager();
		}
		return s_instance;
	}
}
