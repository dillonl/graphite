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
		Reference(Region::SharedPtr region);
		~Reference();

	private:

		std::string m_fasta_path;
	};

}

#endif //GWIZ_REFERENCE_H
