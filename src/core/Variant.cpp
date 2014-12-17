#include "Variant.h"

namespace gwiz
{

	Variant::Variant(const VARIANT_TYPE variant_type,
					 std::string& sequence,
					 position start_position,
					 position end_position) :
		m_variant_type(variant_type),
		m_sequence(sequence),
		m_start_position(start_position),
		m_end_position(end_position)
	{
	}

	Variant::~Variant()
	{
	}
}
