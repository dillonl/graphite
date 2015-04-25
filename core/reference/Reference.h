#ifndef GWIZ_REFERENCE_H
#define GWIZ_REFERENCE_H

#include <stdint.h>
#include <string>
#include <vector>

#include "IReference.h"

namespace gwiz
{

	class Reference : public IReference
	{
	public:
		Reference(std::string& regionString);
		~Reference();

	private:
	};

}

#endif //GWIZ_REFERENCE_H
