#include <boost/graph/graphviz.hpp>

#include "VariantGraph.h"
#include "SNPNode.h"

namespace gwiz
{
	namespace vg
	{

		VariantGraph::VariantGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr) :
			IGraph(referencePtr, variantListPtr),
			m_graph_ptr(std::make_shared< Graph >()),
			m_next_variant_init(false)
		{
		}

		VariantGraph::~VariantGraph()
		{
		}

		VariantGraph::VariantVertexDescriptor VariantGraph::getVertexAtPosition(position referencePosition)
		{
			throw "UNIMPLEMENTED";
			// size_t startIndex = 0;
			// size_t lastIndex = this->m_reference_vertices.size() - 1;
			// while (startIndex <= lastIndex)
			// {
			// 	size_t midIndex = (startIndex + lastIndex) / 2;
			// 	auto midPosition = (*this->m_graph_ptr)[this->m_reference_vertices[midIndex]]->getPosition();
			// 	if (pos > midPosition)
			// 	{
			// 		startIndex = midIndex + 1;
			// 	}
			// 	else if (pos < midPosition)
			// 	{
			// 		lastIndex = midIndex - 1;
			// 	}
			// 	else
			// 	{
			// 		return this->m_reference_vertices[midIndex];
			// 	}
			// }
			// return this->m_reference_vertices[lastIndex];
		}

		void VariantGraph::constructGraph()
		{
			position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
			size_t referenceOffset = 0;
			size_t referenceSize;
			VariantVertexDescriptor referenceVertex;
			std::vector< VariantVertexDescriptor > altAndRefVertices;
			Variant::SharedPtr variantPtr;
			size_t count = 0;
			while (getNextCompoundVariant(variantPtr))
			{
				referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
				if (referenceSize > 0) // there is no reference to add
				{
					ReferenceNode::SharedPtr referenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, referenceSize);
					referenceVertex = addReference(altAndRefVertices, referenceNode);
					altAndRefVertices.clear();
					altAndRefVertices.push_back(referenceVertex);
				}
				size_t variantReferenceSize;
				altAndRefVertices = addVariantVertices(altAndRefVertices, variantPtr, variantReferenceSize);
				referenceOffset += referenceSize + variantReferenceSize;
			}
			referenceSize = this->m_reference_ptr->getRegion()->getEndPosition() - (startPosition + referenceOffset);
			if (referenceSize > 0)
			{
				ReferenceNode::SharedPtr referenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, referenceSize);
				referenceVertex = addReference(altAndRefVertices, referenceNode);
			}
		}

		VariantGraph::VariantVertexDescriptor VariantGraph::addReference(std::vector< VariantGraph::VariantVertexDescriptor >& altAndRefVertices, ReferenceNode::SharedPtr referenceNode)
		{
			// Graph::vertex_descriptor referenceVertex = boost::add_vertex(referenceNode, *m_graph_ptr);
			VariantVertexDescriptor referenceVertex = addReferenceNode(referenceNode);
			// add previous variant and reference Vertices to the referenceVertex
			for (auto iter = altAndRefVertices.begin(); iter != altAndRefVertices.end(); ++iter)
			{
				boost::add_edge((*iter), referenceVertex, *this->m_graph_ptr);
			}
			return referenceVertex;
		}

		std::vector< VariantGraph::VariantVertexDescriptor > VariantGraph::addVariantVertices(std::vector< VariantGraph::VariantVertexDescriptor > altAndRefVertices, Variant::SharedPtr variantPtr, size_t& variantReferenceSize)
		{
			std::vector< VariantVertexDescriptor > vertices;
			for (uint32_t i = 0; i < variantPtr->getAlt().size(); ++i)
			{
				INode::SharedPtr variantNode = IVariantNode::BuildVariantNodes(variantPtr, i);
				vertices.push_back(boost::add_vertex(variantNode, *this->m_graph_ptr));
			}
			size_t referenceOffset = variantPtr->getPosition() - this->m_reference_ptr->getRegion()->getStartPosition();
			ReferenceNode::SharedPtr variantReferenceNode = std::make_shared< ReferenceNode >(this->m_reference_ptr, referenceOffset, variantPtr->getRef()[0].size());
			vertices.push_back(addReferenceNode(variantReferenceNode));
			variantReferenceSize = variantPtr->getRef()[0].size();

			for (auto parentVertexIter = altAndRefVertices.begin(); parentVertexIter != altAndRefVertices.end(); ++parentVertexIter)
			{
				for (auto childVertexIter = vertices.begin(); childVertexIter != vertices.end(); ++childVertexIter)
				{
					boost::add_edge((*parentVertexIter), (*childVertexIter), *this->m_graph_ptr);
				}
			}

			return vertices;
		}

		bool VariantGraph::getNextCompoundVariant(Variant::SharedPtr& variant)
		{
			// the first time we call this function we need to get the first variant
			if (!m_next_variant_init)
			{
				m_next_variant_init = true;
				m_variant_list_ptr->getNextVariant(this->m_next_variant);
			}
			variant = this->m_next_variant; // set the variant to be returned
			if (this->m_next_variant == NULL) { return false; } // if the variant is NULL then return false, this indicates we are past the region of interest

			std::string referenceString = variant->getRef()[0];
			std::vector< Variant::SharedPtr > variants;
			Variant::SharedPtr nextVariant;
			bool variantAdded = false;
			position startPosition = variant->getPosition();

			// loop through variants until the variants stop overlapping.
			// As we loop through build a concatenated reference string
			// that represents the entire overlapped variants reference.
			// for those overlapped variants add them to a vector so
			// a "compound variant" can be generated.
			while(m_variant_list_ptr->getNextVariant(nextVariant))
			{
				position variantEndPosition = (variant->getPosition() + referenceString.size() - 1); // subtract 1 because we are counting starting with the position we are on
				if (variantEndPosition < nextVariant->getPosition())
				{
					break;
				}
				// this is a minor efficiency, even though this is a bit ugly
				// it is more efficient to add the variant after checking that this
				// is a compound variant because it is so rare
				if (!variantAdded)
				{
					variants.push_back(variant); // we will build a compound variant with all these variants
					variantAdded = true;
				}
				std::string nextReferenceString = nextVariant->getRef()[0];
				position nextVariantEndPosition = (nextVariant->getPosition() + nextReferenceString.size());
				// if the next variant has reference at a further position then add it to the end of the referenceString
				if (nextVariantEndPosition > variantEndPosition)
				{
					position referenceDelta = (nextVariantEndPosition - variantEndPosition);
					referenceString += nextReferenceString.substr(referenceDelta);
				}
				variants.push_back(nextVariant); // we will build a compound variant with all these variants
				if (nextVariant->getPosition() < startPosition)
				{
					startPosition = nextVariant->getPosition();
				}
			}
			if (!variants.empty())
			{
				variant = buildCompoundVariant(startPosition, referenceString, variants);
			}
			this->m_next_variant = nextVariant; // set the next variant
			return true;
		}

		Variant::SharedPtr VariantGraph::buildCompoundVariant(const position startPosition, const std::string& referenceString, const std::vector< Variant::SharedPtr >& variants)
		{
			position referenceEndPosition = referenceString.size() + startPosition;
			std::string chrom = variants[0]->getChrom();
			position pos = variants[0]->getPosition();
			std::string id = ".";
			std::string line = chrom + "\t" + std::to_string(pos) + "\t" + id + "\t" + referenceString + "\t";
			// loop over all the variants
			for (auto variantIter = variants.begin(); variantIter != variants.end(); ++variantIter)
			{
				// loop over all the alts in the variants
				for (uint32_t i = 0; i < (*variantIter)->getAlt().size(); ++i)
				{
					std::string altString = (*variantIter)->getAlt()[i];
					// basically we are replacing the variant's reference with the alt within the aggrigated reference (referenceString)
					// it's complicated to explain in words but if you follow the code it isn't too bad
					std::string variantString = referenceString;
					variantString.erase((*variantIter)->getPosition() - startPosition, (*variantIter)->getRef()[0].size());
					variantString.insert((*variantIter)->getPosition() - startPosition, altString);
					line += variantString + ",";
				}
			}
			line.replace(line.size() - 1, 1, "\t\n"); // replace the past comma with a tab
			auto variant = Variant::BuildVariant(line.c_str(), this->m_vcf_parser);
			return variant;
		}

		void VariantGraph::printGraph(const char* path)
		{
			std::ofstream ofs(path);
			boost::write_graphviz(ofs, *this->m_graph_ptr,OurVertexPropertyWriter(*this->m_graph_ptr));
		}

	}// end namespace vg

}// end namespace gwiz
