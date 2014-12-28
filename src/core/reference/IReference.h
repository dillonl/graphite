#ifndef GWIZ_IREFERENCE_H
#define GWIZ_IREFERENCE_H

#include <stdint.h>
#include <memory>

#include "utils/NonCopyable.h"
#include "utils/Types.h"

namespace gwiz
{

	class IReference : private noncopyable
	{
	public:
		typedef std::shared_ptr<IReference> SharedPtr;

		IReference() {}
		virtual ~IReference() {}

	private:
		IReference(const IReference&) = delete;
		IReference& operator=(const IReference&) = delete;
	};
}

#endif // GWIZ_IREFERENCE_H
