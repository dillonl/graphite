#ifndef GRAPHITE_TYPES_H
#define GRAPHITE_TYPES_H

#include <stdint.h>
#include <limits>
#include <string>
#include <vector>

namespace graphite
{
	typedef uint32_t position;
	static position MAX_POSITION = std::numeric_limits< position >::max();
	enum class AlleleCountType { NONE = 0,  NinteyFivePercent = 1, NinteyPercent = 2, EightyPercent = 3, SeventyPercent = 4, AmbiguousPercent = 5 };
	static const std::vector< AlleleCountType > AllAlleleCountTypes = { AlleleCountType::NONE,  AlleleCountType::NinteyFivePercent, AlleleCountType::NinteyPercent, AlleleCountType::EightyPercent, AlleleCountType::SeventyPercent, AlleleCountType::AmbiguousPercent };
	struct AlleleCountTypeHash
	{
		    template <typename T>
			std::size_t operator()(T t) const
			{
				return static_cast<std::size_t>(t);
			}
	};
	static std::string AlleleCountTypeToString(AlleleCountType alleleCountType)
	{
		switch (alleleCountType)
		{
		case AlleleCountType::NinteyFivePercent:
			return "NinteyFivePercent";
		case AlleleCountType::NinteyPercent:
			return "NinteyPercent";
		case AlleleCountType::EightyPercent:
			return "EightyPercent";
		case AlleleCountType::SeventyPercent:
			return "SeventyPercent";
		case AlleleCountType::AmbiguousPercent:
			return "AmbiguousPercent";
		}
		return "NONE";
	}

	static AlleleCountType scoreToAlleleCountType(uint32_t score)
	{
		if (score >= 95)
		{
			return AlleleCountType::NinteyFivePercent;
		}
		else if (score >= 90)
		{
			return AlleleCountType::NinteyPercent;
		}
		else if (score >= 80)
		{
			return AlleleCountType::EightyPercent;
		}
		else if (score >= 70)
		{
			return AlleleCountType::SeventyPercent;
		}
		return AlleleCountType::NONE;
	}
}

#endif //GRAPHITE_TYPES_H
