#ifndef GRAPHITE_IGRAPH_H
#define GRAPHITE_IGRAPH_H

#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/variant/IVariant.h"
#include "core/graph/INode.h"

#include <boost/noncopyable.hpp>

#include <list>
#include <string>
#include <memory>

namespace graphite
{
	class IGraph : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr<IGraph> SharedPtr;

	    IGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
		    m_reference_ptr(referencePtr), m_variant_list_ptr(variantListPtr)
		{
		}
		virtual ~IGraph() {}

	protected:

		virtual void constructGraph() = 0;

		IReference::SharedPtr m_reference_ptr;
		IVariantList::SharedPtr m_variant_list_ptr;
	private:


	};
} // end namespace graphite

#endif // GRAPHITE_IGRAPH_H
