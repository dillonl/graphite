#ifndef GWIZ_IGRAPH_H
#define GWIZ_IGRAPH_H

#include "core/reference/IReference.h"
#include "core/variants/IVariantList.h"
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

	    IGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
		    m_reference_ptr(referencePtr), m_variant_list_ptr(variantListPtr)
		{
		}
		virtual ~IGraph() {}

	protected:
		IGraph() {} // used in tests

		virtual void constructGraph() = 0;
		IReference::SharedPtr m_reference_ptr;
		IVariantList::SharedPtr m_variant_list_ptr;

	private:


	};
} // end namespace gwiz

#endif // GWIZ_IGRAPH_H
