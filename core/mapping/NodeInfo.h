#ifndef GRAPHITE_NODEINFO_H
#define GRAPHITE_NODEINFO_H

namespace graphite
{
    /*
     * NodeInfo stores the length, and variantType of a node. This class is used to generate bed entries based on the graphPathHeaders.
     */
    class NodeInfo
    {
    public:
		typedef std::shared_ptr< NodeInfo > SharedPtr;

        enum VariantType { REF, ALT };

        NodeInfo (int32_t length, VariantType variantType) :
            m_length(length),
            m_variant_type(variantType)
        {
        }

        ~NodeInfo ()
        {
        }

        int32_t getLength () { return m_length; }
        VariantType getVariantType () { return m_variant_type; }

    private:
        int32_t m_length;
        VariantType m_variant_type;
    };
}

#endif
