#pragma once

#include <unordered_map>
#include <vector>
#include <optional>
#include <cmath>
#include <limits>

/**
 * @brief Mappings are interfaces which generalize arrays and hashmaps. This file implements Wrapper for them.
 *
 */
namespace mutil
{
    // ViewIndexProxy helps setting and getting values with the use of indexing
    template<typename V>
    struct ViewIndexProxy {
        V& view;
        typename V::Input index;

        operator typename V::Output()
        {
            return view.get(index);
        }

        ViewIndexProxy& operator=(typename V::Output val){
            view.set(index, val);
            return *this;
        }
    };

    template <typename V>
    struct ConstViewIndexProxy
    {
        const V& view;
        typename V::Input index;

        operator typename V::Output() const
        {
            return view.get(index);
        }
    };

    class Borrowed
    {};

    class Owned
    {};

    /**
     * @brief ArrayMapping is only possible for integer variables as id
     *
     */
    template <class I, class O>
    class ArrayMappingView
    {
    public:
        using Input = I;
        using Output = O;

        using Ownership = Borrowed;

        ArrayMappingView()
            : ArrayMappingView(nullptr, 0)
        {
        }

        ArrayMappingView(O* arr, std::size_t size)
            : ArrayMappingView(arr, size, std::numeric_limits<O>::max())
        {}

        ArrayMappingView(O* arr, std::size_t size, O defaultValue)
            : m_arr(arr)
            , m_size(size)
            , m_defaultValue(defaultValue)
        {}

        void clear(){
            // TODO: might want to use memset
            for (std::size_t i = 0; i < m_size; ++i)
            {
                m_arr[i] = m_defaultValue;
            }
        }

        void clear(Output defaultValue){
            // TODO: might want to use memset
            for (std::size_t i = 0; i < m_size; ++i)
            {
                m_arr[i] = defaultValue;
            }
            m_defaultValue = defaultValue;
        }

        void set(I id, O element)
        {
            m_arr[id] = element;
        }

        O get(I id) const
        {
            return m_arr[id];
        }

        std::optional<Output>
        get_optional(Input id) const
        {
            if (m_arr[id] == m_defaultValue)
            {
                return {};
            }
            return m_arr[id];
        }

        ViewIndexProxy<ArrayMappingView<I,O>>
        operator[](I id) const
        {
            return ViewIndexProxy<ArrayMappingView<I,O>>{*this, id};
        }

        O
        get_unsafe(I id) const
        {
            return m_arr[id];
        }

        void setDefault(O defaultValue){
            if(defaultValue == m_defaultValue){
                return;
            }
            for (std::size_t i = 0; i < m_size; ++i)
            {
                if(m_arr[i] == m_defaultValue){
                    m_arr[i] = defaultValue;
                }
            }
            m_defaultValue = defaultValue;
        }

        O
        getDefault() const
        {
            return m_defaultValue;
        }

        bool
        contains(Input id)
        {
            return m_arr[id] != m_defaultValue;
        }

        O* inner(){
            return m_arr;
        }

    protected:
        O* m_arr;
        std::size_t m_size;
        O m_defaultValue;
    };

    template <class I, class O>
    class VectorMappingView
    {
    public:
        using Input = I;
        using Output = O;

        using Ownership = Owned;

        VectorMappingView()
            : m_vec()
        {}

        VectorMappingView(std::size_t size)
            : m_vec()
        {
            m_vec.resize(size);
        }

        VectorMappingView(std::size_t size, O defaultValue)
            : m_vec()
            , m_defaultValue(defaultValue)
        {
            m_vec.reserve(size);
            for (std::size_t i = 0; i < size; ++i)
            {
                m_vec[i] = defaultValue;
            }
        }

        void
        clear()
        {
            // TODO: might want to use memset
            for (std::size_t i = 0; i < m_vec.size(); ++i)
            {
                m_vec[i] = m_defaultValue;
            }
        }

        void
        clear(Output defaultValue)
        {
            // TODO: might want to use memset
            for (std::size_t i = 0; i < m_vec.size(); ++i)
            {
                m_vec[i] = defaultValue;
            }
            m_defaultValue = defaultValue;
        }

        void
        set(I id, O element)
        {
            m_vec[id] = element;
        }

        O
        get(I id) const
        {
            return m_vec[id];
        }

        std::optional<Output>
        get_optional(Input id) const
        {
            if (m_vec[id] == m_defaultValue)
            {
                return {};
            }
            return m_vec[id];
        }

        // ViewIndexProxy<ArrayMappingView<I, O>>
        // operator[](I id) const
        // {
        //     return ViewIndexProxy<ArrayMappingView<I, O>> { *this, id };
        // }

        O
        get_unsafe(I id) const
        {
            return m_vec[id];
        }

        void
        setDefault(O defaultValue)
        {
            if (defaultValue == m_defaultValue)
            {
                return;
            }
            for (std::size_t i = 0; i < m_vec.size(); ++i)
            {
                if (m_vec[i] == m_defaultValue)
                {
                    m_vec[i] = defaultValue;
                }
            }
            m_defaultValue = defaultValue;
        }

        O
        getDefault() const
        {
            return m_defaultValue;
        }

        bool
        contains(Input id)
        {
            return m_vec[id] != m_defaultValue;
        }

        std::vector<O>&
        inner()
        {
            return m_vec;
        }

    protected:
        std::vector<O> m_vec;
        O m_defaultValue;
    };

    template <class HashMap>
    class HashMapMappingView
    {
    public:
        using Input = typename HashMap::key_type;
        using Output = typename HashMap::mapped_type;

        using Ownership = Owned;

        HashMapMappingView()
            : m_map()
            , m_defaultValue()
        {
        }

        HashMapMappingView(Output defaultValue)
            : m_map()
            , m_defaultValue(defaultValue)
        {
        }

        void clear(){
            m_map.clear();
        }

        void clear(Output defaultValue){
            m_map.clear();
            m_defaultValue = defaultValue;
        }

        void
        set(Input id, Output element)
        {
            m_map[id] = element;
        }

        Output
        get(Input id) const
        {
            auto found = m_map.find(id);
            if (found == m_map.end())
            {
                return m_defaultValue;
            }
            return found->second;
        }

        std::optional<Output>
        get_optional(Input id) const
        {
            auto found = m_map.find(id);
            if (found == m_map.end())
            {
                return {};
            }
            return found->second;
        }

        ViewIndexProxy<HashMapMappingView<HashMap>>
        operator[](Input id)
        {
            return ViewIndexProxy<HashMapMappingView<HashMap>>{*this, id};
        }

        ConstViewIndexProxy<HashMapMappingView<HashMap>>
        operator[](Input id) const
        {
            return ConstViewIndexProxy<HashMapMappingView<HashMap>> { *this, id };
        }

        Output
        get_unsafe(Input id) const
        {
            return m_map.at(id);
        }

        Output
        getDefault() const
        {
            return m_defaultValue;
        }

        bool
        contains(Input id)
        {
            return m_map.find(id) != m_map.end();
        }

        const HashMap& inner() const {
            return m_map;
        }

    protected:
        HashMap m_map;
        Output m_defaultValue;
    };

} // namespace mutil
