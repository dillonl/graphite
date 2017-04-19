#ifndef GRAPHITE_IREFERENCE_H
#define GRAPHITE_IREFERENCE_H

#include <iostream>
#include <stdint.h>
#include <memory>

#include "core/region/Region.h"
#include "core/util/Types.h"

namespace graphite
{

	class IReference : private Noncopyable
	{
	public:
		typedef std::shared_ptr<IReference> SharedPtr;

	    IReference(){}
		virtual ~IReference() {}

		virtual const char* getSequence() { m_sequence.c_str(); }
		virtual size_t getSequenceSize(){ m_sequence.size(); }
		virtual std::string getSequenceFromRegion(Region::SharedPtr regionPtr)
		{
			uint64_t sequenceLength = regionPtr->getEndPosition() - regionPtr->getStartPosition();
			uint32_t startPosition = regionPtr->getStartPosition() - this->m_region->getStartPosition();
			if (regionPtr->getBased() == Region::BASED::ONE)
			{
				startPosition -= 1;;
				sequenceLength += 1;
			}
			std::string sequence = std::string(this->m_sequence.c_str() + startPosition, sequenceLength);
			return sequence;
		}

		Region::SharedPtr getRegion() { return this->m_region; }

	protected:
		Region::SharedPtr m_region;
		std::string m_sequence;

	};
}

#endif // GRAPHITE_IREFERENCE_H
