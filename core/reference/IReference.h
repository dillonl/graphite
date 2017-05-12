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
            //std::cout << "Entered fxn getSequenceFromRegion" << std::endl;
			uint64_t sequenceLength = regionPtr->getEndPosition() - regionPtr->getStartPosition();
            //std::cout << "sequenceLength: " << sequenceLength << std::endl;
			uint32_t startPosition = regionPtr->getStartPosition() - this->m_region->getStartPosition();
            //std::cout << "region start position: " << regionPtr->getStartPosition() << std::endl;
            //std::cout << " reference start position: " << m_region->getStartPosition() << std::endl;
            //std::cout << "startPosition: " << startPosition << std::endl;
			if (regionPtr->getBased() == Region::BASED::ONE)
			{
                std::cout << "Entered if statement BASED::ONE" << std::endl;
				startPosition -= 1;
                std::cout << "startPosition: " << startPosition << std::endl;
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
