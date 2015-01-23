#ifndef GWIZ_IGRAPH_H
#define GWIZ_IGRAPH_H

#include "core/reference/IReference.h"
#include "core/variants/IVariantReader.h"
#include "core/variants/IVariant.h"
#include "core/utils/NonCopyable.h"

#include <list>
#include <string>
#include <memory>

namespace gwiz
{
	enum class GRAPH_TYPE {SPARSE, DENSE};

	class IGraph : private noncopyable
	{
	public:
		typedef std::shared_ptr<IGraph> SharedPtr;

	    IGraph(IReference::SharedPtr referencePtr, IVariantReader::SharedPtr variantReaderPtr) :
		    m_reference_ptr(referencePtr), m_variant_reader_ptr(variantReaderPtr) {}
		virtual ~IGraph() {}

	protected:

		virtual void constructGraph() = 0;
		IReference::SharedPtr m_reference_ptr;
		IVariantReader::SharedPtr m_variant_reader_ptr;

	private:


	};
} // end namespace gwiz

#endif // GWIZ_IGRAPH_H
