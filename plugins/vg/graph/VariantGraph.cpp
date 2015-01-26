#include <boost/graph/graphviz.hpp>

#include "VariantGraph.h"

#include "ReferenceNode.h"
#include "SNPNode.h"

namespace gwiz
{
	namespace vg
	{

		VariantGraph::VariantGraph(IReference::SharedPtr referencePtr, IVariantReader::SharedPtr variantReaderPtr) :
			IGraph(referencePtr, variantReaderPtr),
			m_graph_ptr(std::make_shared< Graph >())
		{
			constructGraph();
		}

		VariantGraph::~VariantGraph()
		{
		}

		void VariantGraph::constructGraph()
		{
			std::vector< std::string > name;
			// boost::dynamic_properties dp;
			// dp.property("label", boost::get(&INode::test, *this->m_graph_ptr));
			position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();

			auto refNode = std::make_shared< ReferenceNode >(this->m_reference_ptr->getSequence(), startPosition, this->m_reference_ptr->getSequenceSize());
			auto v0 = boost::add_vertex(refNode, *m_graph_ptr);

			uint32_t count = 0;
			Variant::SharedPtr variantPtr;
			while (m_variant_reader_ptr->getNextVariant(variantPtr))
			{
				position variantDistance = variantPtr->getPosition() - startPosition;
				refNode->updateLength(variantDistance);
				refNode = std::make_shared< ReferenceNode >(variantPtr->getRef()[0].c_str(), startPosition + variantDistance, variantPtr->getRef()[0].size());
				auto v1 = boost::add_vertex(refNode, *m_graph_ptr);
				boost::add_edge(v0, v1, *m_graph_ptr);
				name.push_back(refNode->sequence);

				for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
				{
					INode::SharedPtr variantNode = std::make_shared< SNPNode >(variantPtr->getAlt()[i].c_str(), variantPtr->getPosition(), variantPtr->getAlt()[i].size());
					auto v2 = boost::add_vertex(variantNode, *m_graph_ptr);
					boost::add_edge(v0, v2, *m_graph_ptr);
					name.push_back(variantNode->sequence);

				}
				++count;
			}
			std::cout << "variant count: " << count << std::endl;

			std::ofstream ofs("out.dot");
			boost::write_graphviz(ofs, *this->m_graph_ptr, boost::make_label_writer(&name[0]));
			// boost::write_graphviz(ofs, *this->m_graph_ptr, boost::make_label_writer(boost::get(&INode::SharedPtr::sequence, *this->m_graph_ptr)));
		}

	}// end namespace vg

}// end namespace gwiz
