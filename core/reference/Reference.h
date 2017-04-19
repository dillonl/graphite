#ifndef GRAPHITE_REFERENCE_H
#define GRAPHITE_REFERENCE_H

#include <stdint.h>
#include <string>
#include <vector>

#include "IReference.h"

namespace graphite
{

	class Reference : public IReference
	{
	public:
		Reference(const std::string& regionSequence, Region::SharedPtr regionPtr);
		~Reference();

	private:
	};

}

#endif //GRAPHITE_REFERENCE_H
