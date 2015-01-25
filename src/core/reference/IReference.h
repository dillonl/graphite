#ifndef GWIZ_IREFERENCE_H
#define GWIZ_IREFERENCE_H

#include <stdint.h>
#include <memory>

#include "core/utils/NonCopyable.h"
#include "core/utils/Types.h"

namespace gwiz
{

	class IReference : private noncopyable
	{
	public:
		typedef std::shared_ptr<IReference> SharedPtr;

		IReference() {}
		virtual ~IReference() {}

		virtual const char* getSequence() = 0;
	private:
		IReference(const IReference&) = delete;
		IReference& operator=(const IReference&) = delete;
	};
}

#endif // GWIZ_IREFERENCE_H
