#include "HTSLibAlignment.h"

namespace graphite
{
	HTSLibAlignment::HTSLibAlignment()
	{
	}

	HTSLibAlignment::~HTSLibAlignment()
	{
	}

	const char* HTSLibAlignment::getSequence()
	{
		std::lock_guard< std::mutex > l(m_sequence_mutex);
		return m_sequence_string.c_str();
	}

	const size_t HTSLibAlignment::getLength()
	{
		std::lock_guard< std::mutex > l(m_sequence_mutex);
		return m_sequence_string.size();
	}

	const void HTSLibAlignment::setSequence(char* seq, uint32_t len)
	{
		std::lock_guard< std::mutex > l(m_sequence_mutex);
		m_sequence_string = std::string(seq, len);
	}

	const void HTSLibAlignment::removeSequence()
	{
		if (true) { return; }
		std::lock_guard< std::mutex > l(m_sequence_mutex);
		--m_sequence_counter;
		if (m_sequence_counter <= 0)
		{
			m_sequence_string.clear();
			m_sequence_counter = 0;
		}
	}

	const void HTSLibAlignment::incrementReferenceCount()
	{
		std::lock_guard< std::mutex > l(m_sequence_mutex);
		++m_sequence_counter;
	}
}
