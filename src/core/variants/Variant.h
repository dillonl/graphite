#ifndef GWIZ_VARIANT_H
#define GWIZ_VARIANT_h

#include <string>
#include <vector>

#include "reference/IReference.h"
#include "IVariant.h"

namespace gwiz
{

	enum class VARIANT_TYPE
	{
		SNP,
		INS,
		DEL,
		DUP,
		INV,
	};

	class Variant : public IVariant
	{
	public:
		Variant(const VARIANT_TYPE variant_type,
				std::string& sequence,
				position start_position,
			    position end_position);
		~Variant();

		VARIANT_TYPE getVariantType() { return m_variant_type; }
		position getStartPosition() { return m_start_position; }
		position getEndPosition() { return m_end_position; }

	private:
		VARIANT_TYPE m_variant_type;
		std::string m_sequence;
		position m_start_position;
		position m_end_position;
	};

}

#endif //GWIZ_VARIANT_h
