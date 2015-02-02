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
			m_graph_ptr(std::make_shared< Graph >()),
			m_next_variant_init(false)
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
				for (auto varAltIter = (*variantIter)->getAlt().begin(); varAltIter != (*variantIter)->getAlt().end(); ++varAltIter)
				{
					// basically we are replacing the variant's reference with the alt within the aggrigated reference (referenceString)
					// it's complicated to explain in words but if you follow the code it isn't too bad
					std::string altString = (*varAltIter);
					std::string variantString = referenceString;
					variantString.erase((*variantIter)->getPosition() - startPosition, (*variantIter)->getRef()[0].size());
					variantString.insert((*variantIter)->getPosition() - startPosition, altString);

					line += variantString + ",";
				}
			}

			line.replace(line.size() - 1, 1, "\t"); // replace the past comma with a tab
			return Variant::BuildVariant(line.c_str(), this->m_vcf_parser);
		}

		void VariantGraph::printGraph(const char* path)
		{
			std::ofstream ofs(path);
			boost::write_graphviz(ofs, *this->m_graph_ptr,OurVertexPropertyWriter(*this->m_graph_ptr));
		}

	}// end namespace vg

}// end namespace gwiz
