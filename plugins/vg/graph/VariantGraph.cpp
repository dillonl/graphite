#include <boost/graph/graphviz.hpp>

#include "VariantGraph.h"

#include "ReferenceNode.h"
#include "SNPNode.h"

namespace gwiz
{
	namespace vg
	{

		VariantGraph::VariantGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
			IGraph(referencePtr, variantListPtr),
			m_graph_ptr(std::make_shared< Graph >())
		{
			constructGraph();
		}

		VariantGraph::~VariantGraph()
		{
		}

		void VariantGraph::constructGraph()
		{
			// std::lock_guard< std::mutex > lock(this->m_graph_mutex); // m_graph_mutex lock will release when it falls out of scope
			position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
			size_t referenceOffset = 0;
			Variant::SharedPtr variantPtr;
			std::vector< Graph::vertex_descriptor > altAndRefVertices;
			size_t count = 0;
			Graph::vertex_descriptor referenceVertex;
			size_t referenceSize;
			while (m_variant_list_ptr->getNextVariant(variantPtr))
			{
				referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
				if (referenceSize > 0)
				{
					std::cout << "ref > 0: " << variantPtr->getPosition() << std::endl;
					std::cout << "refsize: " << referenceSize << std::endl;
					auto referenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, referenceSize);
					referenceVertex = boost::add_vertex(referenceNode, *m_graph_ptr);
					// add previous variant and reference Vertices to the referenceVertex
					for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
					{
						boost::add_edge((*iter), referenceVertex, *this->m_graph_ptr);
					}

					altAndRefVertices.clear(); // clear out the alt and ref vertices
					// get next variants and add edges from the ref to the variants
					for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
					{
						INode::SharedPtr variantNode = IVariantNode::BuildVariantNodes(variantPtr, i);
						auto variantVertex = boost::add_vertex(variantNode, *this->m_graph_ptr);
						boost::add_edge(referenceVertex, variantVertex, *this->m_graph_ptr);
						altAndRefVertices.push_back(variantVertex); // add this vertex so the next reference can add an edge between itself and this variant
					}

					//adding the reference node from the variantPtr
					ReferenceNode::SharedPtr variantReferenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset + referenceSize, variantPtr->getRef()[0].size());
					auto variantReferenceVertex = boost::add_vertex(variantReferenceNode, *m_graph_ptr);
					boost::add_edge(referenceVertex, variantReferenceVertex, *this->m_graph_ptr);
					altAndRefVertices.push_back(variantReferenceVertex);
				}
				else
				{
					std::cout << "ref == 0: " << variantPtr->getPosition() << std::endl;
					std::vector< Graph::vertex_descriptor > tmpAltAndRefVertices;
					for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
					{
						INode::SharedPtr variantNode = IVariantNode::BuildVariantNodes(variantPtr, i);
						auto variantVertex = boost::add_vertex(variantNode, *this->m_graph_ptr);
						for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
						{
							boost::add_edge((*iter), variantVertex, *this->m_graph_ptr);
						}
						tmpAltAndRefVertices.push_back(variantVertex); // add this vertex so the next reference can add an edge between itself and this variant
					}
					ReferenceNode::SharedPtr variantReferenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset + referenceSize, variantPtr->getRef()[0].size());
					auto variantReferenceVertex = boost::add_vertex(variantReferenceNode, *m_graph_ptr);
					for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
					{
						boost::add_edge((*iter), variantReferenceVertex, *this->m_graph_ptr);
					}
					tmpAltAndRefVertices.push_back(variantReferenceVertex);
					referenceOffset += referenceSize + variantPtr->getRef()[0].size();
					altAndRefVertices = tmpAltAndRefVertices;
				}
				referenceOffset += referenceSize + variantPtr->getRef()[0].size();
			}

			std::cout << "almost done" << std::endl;
			// if the referenceSize == 0
			uint32_t endPosition = (startPosition + referenceOffset);
			if (this->m_reference_ptr->getRegion()->getEndPosition() > endPosition)
			{
				referenceSize = this->m_reference_ptr->getRegion()->getEndPosition() - endPosition;
				auto referenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, referenceSize);
				referenceVertex = boost::add_vertex(referenceNode, *m_graph_ptr);
				// add previous variant and reference Vertices to the referenceVertex
				for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
				{
					boost::add_edge((*iter), referenceVertex, *this->m_graph_ptr);
				}
			}
		}

		bool getCompoundNode(Variant::SharedPtr variant)
		{
			std::vector< Variant::SharedPtr > variants;
			Variant::SharedPtr nextVariant;
			/*
			while(m_variant_list_ptr->getNextVariant(nextVariant) && doVariantsOverlap())
			{

			}

			if (m_variant_list_ptr->getNextVariant(nextVariant))
			{

			}
			*/
			return false;
		}

		Variant::SharedPtr buildCompoundNode(std::vector< Variant::SharedPtr > variants)
		{
			return std::make_shared< Variant >();
		}

		void VariantGraph::printGraph(const char* path)
		{
			std::ofstream ofs(path);
			boost::write_graphviz(ofs, *this->m_graph_ptr,OurVertexPropertyWriter(*this->m_graph_ptr));
		}

		/*
		void VariantGraph::constructGraph()
		{
			position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
			size_t referenceOffset = 0;
			Variant::SharedPtr variantPtr;
			m_variant_reader_ptr->getNextVariant(variantPtr);
			ReferenceNode::SharedPtr refNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, variantPtr->getPosition() - startPosition);
			auto referenceVertex = boost::add_vertex(refNode, *m_graph_ptr);

			// -- Should be converted into a do while loop

			while (m_variant_reader_ptr->getNextVariant(variantPtr))
			{
				 // ReferenceNode::SharedPtr refNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, )
				for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
				{
					// IVariantNode::SharedPtr variantNode = IVariantNode::BuildVariantNodes(variantPtr, i);
					auto variantVertex = boost::add_vertex(IVariantNode::BuildVariantNodes(variantPtr, i), *this->m_graph_ptr);
					boost::add_edge(referenceVertex, variantVertex, *this->m_graph_ptr);
				}
			}
		}
		*/

		/*
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
				std::cout << refNode->sequence << std::endl;
				exit(0);

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
		*/
	}// end namespace vg

}// end namespace gwiz
