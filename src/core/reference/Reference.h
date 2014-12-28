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
		Reference(const std::string& fasta_path);
		~Reference();

	private:

		std::string m_fasta_path;
	};

}

#endif //GWIZ_REFERENCE_H
