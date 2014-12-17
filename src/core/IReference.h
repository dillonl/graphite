#ifndef GWIZ_IREFERENCE_H
#define GWIZ_IREFERENCE_H

#include <stdint.h>

#include "utils/NonCopyable.h"

namespace gwiz
{
	typedef uint32_t position;

	class IReference : private noncopyable
	{
	public:
		IReference() {}
		virtual ~IReference() {}

	private:
		IReference(const IReference&) = delete;
		IReference& operator=(const IReference&) = delete;
	};
}

#endif // GWIZ_IREFERENCE_H
