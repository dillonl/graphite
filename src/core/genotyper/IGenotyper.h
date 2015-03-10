#ifndef GWIZ_GENOTYPER_H
#define GWIZ_GENOTYPER_H

#include "IGenotyperVariant.h"

#include <memory>
#include <iostream>
#include <fstream>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class IGenotyper : private boost::noncopyable
	{
	public:
		static IGenotyper* Instance()
		{
			static IGenotyper* s_genotyper;
			if (s_genotyper == nullptr)
			{
				s_genotyper = new IGenotyper();
			}
			return s_genotyper;
		}

		IGenotyperVariant::SharedPtr generateVariant(position pos)
		{
			auto genotyperVariant = std::make_shared< IGenotyperVariant >(pos);
			this->m_variant_mutex.lock();
			this->m_genotyper_variant.push_back(genotyperVariant);
			this->m_variant_mutex.unlock();
			return genotyperVariant;
		}

		void printVariants()
		{
this->m_genotyper_variant.sort([](const IGenotyperVariant::SharedPtr& a, const IGenotyperVariant::SharedPtr& b) -> bool {
return a->getPosition() > b->getPosition();
});
			std::ofstream outputFile;
			outputFile.open("graph_output.txt");
			for(auto genotypeVariantPtr : this->m_genotyper_variant)
			{
                std::string referenceString = "";
                std::string variantString = "";
				for (auto allele : genotypeVariantPtr->getGenotyperAlleles())
				{
					if (allele->getType() == GenotyperAllele::Type::REFERENCE)
					{
						referenceString = allele->getSequence() + "<" + std::to_string(allele->getReadCount()) + ">\t";
					}
					else
					{
						variantString += allele->getSequence() + "<" + std::to_string(allele->getReadCount()) + ">\t";
					}
				}
                outputFile << genotypeVariantPtr->getPosition() << "\t" << referenceString << variantString << std::endl;
			}
		}

	protected:
		std::list< IGenotyperVariant::SharedPtr > m_genotyper_variant;
		std::mutex m_variant_mutex;

	private:
	    IGenotyper()
		{
		}

		virtual ~IGenotyper()
		{
		}

		static IGenotyper* s_genotyper;
	};
}

#endif //GWIZ_GENOTYPER_H
