#ifndef GWIZ_TYPES_H
#define GWIZ_TYPES_H

#include <stdint.h>
#include <limits>

namespace gwiz
{
	typedef uint32_t position;
	static position MAX_POSITION = std::numeric_limits< position >::max();
}

#endif //GWIZ_TYPES_H
