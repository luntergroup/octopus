//
//  matrix_map.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef matrix_map_h
#define matrix_map_h

#include <vector>
#include <unordered_map>
#include <functional>
#include <cstddef>
#include <utility>
#include <iterator>
#include <type_traits>

template <typename Key1,
          typename Key2,
          typename T,
          typename Hash1 = std::hash<Key1>,
          typename Hash2 = std::hash<Key2>,
          typename KeyEqual1 = std::equal_to<Key1>,
          typename KeyEqual2 = std::equal_to<Key2>
> class MatrixMap
{
public:
    using size_type = std::size_t;
    
    MatrixMap()  = default;
    ~MatrixMap() = default;
    
    template <typename K>
    MatrixMap(K&& key2s) : key2s_ {std::forward<T>(key2s)}
    {
        key2_indicies_.reserve(key2s_.size());
    }
    
    MatrixMap(const MatrixMap& other)
    :
    key2s_ {other.key2s_},
    values_ {other.values_}
    {
        this->generate_incice_map(other.key2_indicies_);
    }
    
    MatrixMap& operator=(const MatrixMap& other)
    {
        if (&other == this) {
            return *this;
        }
        
        key2s_  = other.key2s_;
        values_ = other.values_;
        
        this->generate_incice_map(other.key2_indicies_);
        
        return *this;
    }
    
    MatrixMap(MatrixMap&& other)
    :
    key2s_ {std::move(other.key2s_)},
    values_ {std::move(other.values_)}
    {
        this->generate_incice_map(other.key2_indicies_);
    }
    
    MatrixMap& operator=(MatrixMap&& other)
    {
        if (&other == this) {
            return *this;
        }
        
        key2s_  = std::move(other.key2s_);
        values_ = std::move(other.values_);
        
        this->generate_incice_map(other.key2_indicies_);
        
        return *this;
    }
    
    T& operator()(const Key1& key1, const Key2& key2)
    {
        return values_.at(key1)[key2_indicies_.at(key2)];
    }
    
    const T& operator()(const Key1& key1, const Key2& key2) const
    {
        return values_.at(key1)[key2_indicies_.at(key2)];
    }
    
    bool empty1() const noexcept
    {
        return values_.empty();
    }
    
    bool empty2() const noexcept
    {
        return key2s_.empty();
    }
    
    size_type size1() const noexcept
    {
        return values_.size();
    }
    
    size_type size2() const noexcept
    {
        return key2s_.size();
    }
    
    void reserve(size_type n1)
    {
        values_.reserve(n1);
    }
    
    void reserve(size_type n1, size_type n2)
    {
        values_.reserve(n1);
        key2s_.reserve(n2);
        key2_indicies_.reserve(n2);
    }
    
    void clear() noexcept
    {
        key2s_.clear();
        values_.clear();
        key2_indicies_.clear();
    }
    
    template <typename InputIt>
    bool assign_key2s(InputIt first, InputIt last)
    {
        key2s_.assign(first, last);
        
        this->regenerate_indice_map();
        
        if (!values_.empty()) {
            values_.clear();
            return true;
        }
        
        return false;;
    }
    
    bool push_back_key2(const Key2& key)
    {
        this->push_back_reallocate(key);
        
        if (!values_.empty()) {
            values_.clear();
            return true;
        }
        
        return false;
    }
    
    template <typename InputIt>
    void insert_at(const Key1& key, InputIt first, InputIt last)
    {
        if (std::distance(first, last) != this->size2()) {
            throw std::out_of_range {"MatrixMap::insert_at called with value range of different"
                " length to Key2 range in this MatrixMap"};
        }
        values_[key].assign(first, last);
    }
    
    template <typename InputIt>
    bool insert_each(const Key2& key, InputIt first, InputIt last)
    {
        if (std::distance(first, last) != this->size1()) {
            throw std::out_of_range {"MatrixMap::insert_each called with value range of different"
                " length to Key1 range in this MatrixMap"};
        }
        
        if (key2_indicies_.count(key) != 0) {
            return false;
        }
        
        for (auto& p : values_) {
            p.second.push_back(*first);
            ++first;
        }
        
        this->push_back_reallocate(key);
        
        return true;
    }
    
private:
    struct Key2RefHash
    {
        std::size_t operator()(std::reference_wrapper<const Key2> k) const
        {
            return Hash2()(k.get());
        }
    };
    
    struct Key2RefEqual
    {
        bool operator()(const Key2& lhs, std::reference_wrapper<const Key2> rhs) const
        {
            return KeyEqual2()(lhs, rhs.get());
        }
    };
    
    using Key2ContainerType = std::vector<Key2>;
    using ValueContainerType = std::vector<T>;
    using ValueMap = std::unordered_map<Key1, ValueContainerType, Hash1, KeyEqual1>;
    using IndexSizeType = typename ValueMap::mapped_type::size_type;
    using IndiceMap = std::unordered_map<std::reference_wrapper<const Key2>,
    IndexSizeType, Key2RefHash, Key2RefEqual>;
    
    using Key2Iterator  = typename Key2ContainerType::const_iterator;
    using ValueIterator = typename ValueContainerType::const_iterator;
    
public:
    
    class Iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = std::pair<const Key2&, const T&>;
        using difference_type   = typename std::iterator_traits<Key2Iterator>::difference_type;
        using reference         = value_type;
        using pointer           = value_type*;
        
        Iterator() = delete;
        
        explicit Iterator(Key2Iterator key2_itr, ValueIterator value_itr)
        : key2_itr_ {key2_itr}, value_itr_ {value_itr} {}
        
        ~Iterator() = default;
        
        Iterator& operator++()
        {
            key2_itr_++;
            value_itr_++;
            return *this;
        }
        
        value_type operator*() const
        {
            return std::make_pair(std::ref(*key2_itr_), std::ref(*value_itr_));
        }
        
        friend bool operator==(const Iterator& lhs, const Iterator& rhs)
        {
            return lhs.key2_itr_ == rhs.key2_itr_;
        }
        
        friend bool operator!=(const Iterator& lhs, const Iterator& rhs)
        {
            return !operator==(lhs, rhs);
        }
    private:
        Key2Iterator key2_itr_;
        ValueIterator value_itr_;
    };
    
private:
    
    Key2ContainerType key2s_;
    ValueMap values_;
    IndiceMap key2_indicies_;
    
    void generate_incice_map(const IndiceMap& other)
    {
        key2_indicies_.reserve(key2s_.size());
        
        for (const auto& key : key2s_) {
            key2_indicies_.emplace(key, other.at(key));
        }
    }
    
    void regenerate_indice_map()
    {
        key2_indicies_.clear();
        
        std::size_t i {0};
        
        for (const auto& key : key2s_) {
            key2_indicies_.emplace(key, i);
        }
    }
    
    void push_back_reallocate(const Key2& key)
    {
        if (key2s_.capacity() == key2s_.size()) {
            key2s_.push_back(key);
            this->regenerate_indice_map();
        } else {
            key2s_.push_back(key);
            key2_indicies_.emplace(key2s_.back(), key2s_.size() - 1);
        }
    }
    
public:
    
    Iterator begin(const Key2& key) const
    {
        return Iterator {std::begin(key2s_), std::begin(values_.at(key))};
    }
    
    Iterator end(const Key2& key) const
    {
        return Iterator {std::end(key2s_), std::end(values_.at(key))};
    }
    
    Iterator cbegin(const Key2& key) const
    {
        return begin(key);
    }
    
    Iterator cend(const Key2& key) const
    {
        return end(key);
    }
    
    class Proxy
    {
    public:
        Proxy() = delete;
        ~Proxy() = default;
        Proxy(Iterator begin, Iterator end) : begin_ {begin}, end_ {end} {}
        Iterator begin() const { return begin_; }
        Iterator end() const { return end_; }
        Iterator cbegin() const { return begin_; }
        Iterator cend() const { return end_; }
    private:
        Iterator begin_, end_;
    };
    
    Proxy operator[](const Key1& key) const
    {
        return Proxy {this->begin(key), this->end(key)};
    }
};

#endif /* matrix_map_h */
