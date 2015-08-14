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
		Reference(std::string& regionString);
		~Reference();

	private:
	};

}

#endif //GRAPHITE_REFERENCE_H
