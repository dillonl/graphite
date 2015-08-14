#ifndef GRAPHITE_TYPES_H
#define GRAPHITE_TYPES_H

#include <stdint.h>
#include <limits>

namespace graphite
{
	typedef uint32_t position;
	static position MAX_POSITION = std::numeric_limits< position >::max();
}

#endif //GRAPHITE_TYPES_H
